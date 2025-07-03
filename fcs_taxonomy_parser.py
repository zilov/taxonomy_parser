#!/usr/bin/env python3

import sys
import argparse
import logging
from collections import defaultdict
from typing import Generator, Tuple, Dict
import pandas as pd
import gzip
import os

def setup_logger():
    """Configure and return a logger for the application."""
    logger = logging.getLogger('taxonomy_parser')
    logger.setLevel(logging.INFO)
    
    # Create console handler
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    
    # Add handler to logger
    logger.addHandler(ch)
    return logger

# Initialize logger
logger = setup_logger()

def parse_taxonomy_file(file_path):
    """Parse the taxonomy file and extract scaffold-taxa relationships with coverage."""
    scaffold_taxa_lengths = defaultdict(lambda: defaultdict(int))
    
    with open(file_path, 'r') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 10:  # Ensure we have enough fields
                continue
            
            # Extract scaffold ID and remove coordinates if present
            full_scaffold_id = fields[0]
            scaffold_id = full_scaffold_id.split('~')[0] if '~' in full_scaffold_id else full_scaffold_id
            
            # Extract tax IDs and their corresponding coverages
            tax_id_positions = [6, 12, 18, 24]
            cvg_positions = [9, 15, 21, 27]
            
            for tax_pos, cvg_pos in zip(tax_id_positions, cvg_positions):
                if tax_pos < len(fields) and cvg_pos < len(fields):
                    tax_id = fields[tax_pos]
                    if tax_id == 'n/a':
                        continue
                    
                    try:
                        coverage = int(fields[cvg_pos])
                        scaffold_taxa_lengths[scaffold_id][tax_id] += coverage
                    except (ValueError, IndexError):
                        continue
    
    return scaffold_taxa_lengths

def fasta_reader_yield(path_to_fasta_file: str) -> Generator[Tuple[str, str], None, None]:
    """Yield tuples of (header, sequence) from a FASTA file (plain text or gzipped)."""
    
    header = None
    seq = []
    
    # Determine if file is gzipped based on extension
    is_gzipped = path_to_fasta_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with opener(path_to_fasta_file, mode) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def get_fasta_lengths(fasta_file: str) -> Dict[str, int]:
    """Create a dictionary mapping sequence headers to their lengths."""
    header_length_dict = {}
    for header, sequence in fasta_reader_yield(fasta_file):
        scaffold_id = header.split()[0]
        header_length_dict[scaffold_id] = len(sequence)
    return header_length_dict

def generate_summary(scaffold_taxa_lengths):
    """Generate summary of taxonomic matches for each scaffold."""
    summary = []
    
    for scaffold_id, taxa_lengths in scaffold_taxa_lengths.items():
        # Sort taxa by length in descending order
        sorted_taxa = sorted(taxa_lengths.items(), key=lambda x: x[1], reverse=True)
        
        # Generate rows with rank
        for rank, (tax_id, length) in enumerate(sorted_taxa, 1):
            summary.append((scaffold_id, tax_id, length, rank))
    
    return summary

def parse_ranked_lineage(lineage_file):
    """Parse the NCBI ranked lineage file and return a DataFrame with lineage information."""
    columns = ['taxid', 'tax_name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
    
    # Fix the separator pattern - use raw string to avoid escape sequence warning
    lineage_df = pd.read_csv(lineage_file, sep=r'\t\|\t', engine='python', names=columns)
    lineage_df[columns[-1]] = lineage_df[columns[-1]].str.rstrip('\t|')
    
    # Clean all string columns
    for col in columns[1:]:  # Skip taxid
        lineage_df[col] = lineage_df[col].astype(str).str.strip()
    
    # Convert taxid to string and set as index for faster lookups
    lineage_df['taxid'] = lineage_df['taxid'].astype(str)
    lineage_df.set_index('taxid', inplace=True)
    
    logger.info(f"Loaded {len(lineage_df)} taxonomic records from lineage database")
    
    return lineage_df

def is_target_taxon(tax_id, lineage_df, target_taxa):
    """Check if a tax ID belongs to the target taxa based on its lineage."""
    try:
        tax_id_str = str(tax_id).strip()
        if tax_id_str not in lineage_df.index or not target_taxa:
            return False
        
        row = lineage_df.loc[tax_id_str]
        
        for level, target_value in target_taxa.items():
            if level in row.index and pd.notna(row[level]):
                actual_value = str(row[level]).lower().strip()
                target_value_clean = target_value
                
                # Skip 'nan' values
                if actual_value == 'nan':
                    continue
                
                logger.debug(f"Comparing {level}: '{actual_value}' vs '{target_value_clean}' for tax_id {tax_id_str}")
                
                if actual_value == target_value_clean:
                    logger.debug(f"Target match found: {tax_id_str} matches {level}={target_value_clean}")
                    return True
        
        return False
    except Exception as e:
        logger.debug(f"Error checking target taxon for {tax_id}: {e}")
        return False

def write_output(summary, output_file, fasta_lengths=None, lineage_df=None, target_taxa=None):
    """Write output to a file with all columns, using None for unavailable data."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("scaffold_id,tax_id,covered_length,top_n,percentage_covered,is_target\n")
        
        for row in summary:
            scaffold_id, tax_id, covered_length, rank = row
            
            # Calculate percentage if fasta lengths are available
            percentage = None
            if fasta_lengths and scaffold_id in fasta_lengths:
                total_length = fasta_lengths[scaffold_id]
                if total_length > 0:
                    percentage = (covered_length / total_length) * 100
            
            # Check if taxon is a target
            is_target = None
            if lineage_df is not None and target_taxa:
                is_target = is_target_taxon(tax_id, lineage_df, target_taxa)
            
            # Format the values
            percentage_str = f"{percentage:.2f}" if percentage is not None else "None"
            is_target_str = str(is_target) if is_target is not None else "None"
            
            # Write the output
            f.write(f"{scaffold_id},{tax_id},{covered_length},{rank},{percentage_str},{is_target_str}\n")

def write_non_target_scaffolds(summary_file, output_file, lineage_df):
    """Write non-target scaffolds with lineage information of top1 match using the summary output file."""
    # Read the summary file to identify non-target scaffolds
    df = pd.read_csv(summary_file)
    
    # Group by scaffold_id and check if any row has is_target=True
    scaffold_groups = df.groupby('scaffold_id')
    non_target_scaffolds = []
    
    for scaffold_id, group in scaffold_groups:
        # Check if this scaffold has any target taxa matches
        has_target = group['is_target'].any()
        
        if not has_target:
            # Get the top1 match (top_n == 1)
            top1_rows = group[group['top_n'] == 1]
            if not top1_rows.empty:
                top1_row = top1_rows.iloc[0]
                non_target_scaffolds.append(top1_row)
    
    # Write the non-target scaffolds file
    with open(output_file, 'w') as f:
        # Write header
        f.write("scaffold_id,top_1_taxa,percentage_covered,species,genus,family,order,class,phylum,kingdom\n")
        
        for _, row in enumerate(non_target_scaffolds):
            scaffold_id = row['scaffold_id']
            tax_id = row['tax_id']
            percentage_covered = row['percentage_covered']
            
            # Get lineage information
            lineage_info = ["None"] * 7  # species, genus, family, order, class, phylum, kingdom
            try:
                tax_id_str = str(tax_id).strip()
                if tax_id_str in lineage_df.index:
                    lineage_row = lineage_df.loc[tax_id_str]
                    lineage_columns = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
                    for i, col in enumerate(lineage_columns):
                        if col in lineage_row.index and pd.notna(lineage_row[col]) and str(lineage_row[col]).strip() != 'nan':
                            lineage_info[i] = str(lineage_row[col]).strip()
            except Exception as e:
                logger.debug(f"Error getting lineage for {tax_id}: {e}")
            
            # Write the output
            lineage_str = ",".join(lineage_info)
            f.write(f"{scaffold_id},{tax_id},{percentage_covered},{lineage_str}\n")

def main(tax_file=None, fasta_file=None, lineage_file=None, target_taxa=None, outdir=None):
    """Main function to run the parser."""
    
    # Initialize variables
    fasta_lengths = {}
    lineage_df = None
    scaffold_taxa_lengths = defaultdict(lambda: defaultdict(int))
    
    # Process fasta file if provided
    if fasta_file:
        logger.info(f"Processing fasta file: {fasta_file}")
        fasta_lengths = get_fasta_lengths(fasta_file)
        if not fasta_lengths:
            logger.error(f"No sequences found in FASTA file {fasta_file}. The file may be empty or in an incorrect format.")
            raise ValueError(f"No sequences found in FASTA file {fasta_file}. The file may be empty or in an incorrect format.")
    
    # Process lineage file if provided
    if lineage_file:
        logger.info(f"Processing lineage database: {lineage_file}")
        lineage_df = parse_ranked_lineage(lineage_file)
    
    if target_taxa:
        logger.info(f"Target taxa filters: {target_taxa}")
            

    logger.info(f"Processing taxonomy file: {tax_file}")
    scaffold_taxa_lengths = parse_taxonomy_file(tax_file)
    
    # Generate summary
    summary = generate_summary(scaffold_taxa_lengths)
    
    # Check for scaffold IDs in taxonomy file but not in FASTA file
    if fasta_file:
        missing_scaffolds = [sid for sid in scaffold_taxa_lengths.keys() if sid not in fasta_lengths]
        if missing_scaffolds:
            logger.warning(f"Scaffold IDs from taxonomy file not found in FASTA file: {missing_scaffolds[:5]}...")
            raise ValueError(f"Scaffold IDs from taxonomy file not found in FASTA file: {missing_scaffolds[:5]}...")
    
    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    
    # Generate output filename
    tax_basename = os.path.basename(tax_file)
    output_filename = f"{tax_basename}.summary.csv"
    output_file = os.path.join(outdir, output_filename)

    logger.info(f"Writing output to: {output_file}")
    
    # Write output
    write_output(summary, output_file, fasta_lengths, lineage_df, target_taxa)
    
    # Write non-target scaffolds if lineage and target taxa are available
    if lineage_df is not None and target_taxa:
        non_target_filename = f"{tax_basename}.non_target.csv"
        non_target_file = os.path.join(outdir, non_target_filename)
        logger.info(f"Writing non-target scaffolds to: {non_target_file}")
        write_non_target_scaffolds(output_file, non_target_file, lineage_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse taxonomy file and extract scaffold-taxa relationships with coverage.")
    parser.add_argument('-t', '--tax_results', help='Path to the taxonomy file to parse', required=True)
    parser.add_argument('-f', '--fasta', help='Path to the fasta file to process')
    parser.add_argument('-r', '--rankedlineage_db', help='Path to NCBI ranked lineage database file')
    parser.add_argument('--target_taxa', nargs='+', help='Target taxa in format taxon:value (e.g., order:coleoptera family:primates)')
    parser.add_argument('--log', help='Path to log file (if not provided, logs to stderr)', default=None)
    parser.add_argument('-o', '--outdir', help='Output directory for results (default: current directory)', default=None, required=True)
    args = parser.parse_args()
    
    # Configure file logging if requested
    if args.log:
        file_handler = logging.FileHandler(args.log)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levellevel)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Parse target taxa
    target_taxa = {}
    if args.target_taxa:
        for item in args.target_taxa:
            try:
                level, value = item.split(':')
                target_taxa[level.lower().strip()] = value.lower().strip()
                logger.info(f"Added target taxa filter: {level.lower().strip()} = {value.lower().strip()}")
            except ValueError:
                logger.warning(f"Invalid target taxa format '{item}'. Expected format is 'taxon:value'")
    
    main(tax_file=args.tax_results, fasta_file=args.fasta, lineage_file=args.rankedlineage_db, target_taxa=target_taxa, outdir=args.outdir)