#!/usr/bin/env python3

import sys
import argparse
import logging
from collections import defaultdict
from typing import Generator, Tuple, Dict
import pandas as pd
import gzip

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
    """Parse the NCBI ranked lineage file and create a mapping of tax IDs to lineage information."""
    columns = ['taxid', 'tax_name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
    
    lineage_df = pd.read_csv(lineage_file, sep='\t\|\t', engine='python', names=columns)
    lineage_df[columns[-1]] = lineage_df[columns[-1]].str.rstrip('\t|')
    
    # Convert DataFrame to dictionary for easy lookups
    lineage_dict = {}
    for _, row in lineage_df.iterrows():
        taxid = row['taxid']
        lineage_dict[taxid] = {col: row[col] for col in columns[1:]}
    
    return lineage_dict

def is_target_taxon(tax_id, lineage_dict, target_taxa):
    """Check if a tax ID belongs to the target taxa based on its lineage."""
    if tax_id not in lineage_dict or not target_taxa:
        return False
    
    lineage = lineage_dict[tax_id]
    for level, target_value in target_taxa.items():
        if level in lineage and lineage[level].lower() == target_value.lower():
            return True
    
    return False

def write_output(summary, fasta_lengths=None, lineage_dict=None, target_taxa=None):
    """Write output with all columns, using None for unavailable data."""
    # Print header
    print("scaffold_id,tax_id,covered_length,top_n,percentage_covered,is_target")
    
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
        if lineage_dict and target_taxa:
            is_target = is_target_taxon(tax_id, lineage_dict, target_taxa)
        
        # Format the values
        percentage_str = f"{percentage:.2f}" if percentage is not None else "None"
        is_target_str = str(is_target) if is_target is not None else "None"
        
        # Print the output
        print(f"{scaffold_id},{tax_id},{covered_length},{rank},{percentage_str},{is_target_str}")

def main(tax_file=None, fasta_file=None, lineage_file=None, target_taxa=None):
    """Main function to run the parser."""
    
    # Initialize variables
    fasta_lengths = {}
    lineage_dict = {}
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
        lineage_dict = parse_ranked_lineage(lineage_file)
    
    if target_taxa:
        logger.info(f"Target taxa filters: {target_taxa}")
    
    # Process taxonomy file if provided
    if tax_file:
        logger.info(f"Processing taxonomy file: {tax_file}")
        scaffold_taxa_lengths = parse_taxonomy_file(tax_file)
    
    # Generate summary
    summary = generate_summary(scaffold_taxa_lengths)
    
    # Check for scaffold IDs in taxonomy file but not in FASTA file
    if tax_file and fasta_file:
        missing_scaffolds = [sid for sid in scaffold_taxa_lengths.keys() if sid not in fasta_lengths]
        if missing_scaffolds:
            logger.warning(f"Scaffold IDs from taxonomy file not found in FASTA file: {missing_scaffolds[:5]}...")
            raise ValueError(f"Scaffold IDs from taxonomy file not found in FASTA file: {missing_scaffolds[:5]}...")
    
    # Write output
    write_output(summary, fasta_lengths, lineage_dict, target_taxa)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse taxonomy file and extract scaffold-taxa relationships with coverage.")
    parser.add_argument('-t', '--tax_results', help='Path to the taxonomy file to parse', required=True)
    parser.add_argument('-f', '--fasta', help='Path to the fasta file to process')
    parser.add_argument('-r', '--rankedlineage_db', help='Path to NCBI ranked lineage database file')
    parser.add_argument('--target_taxa', nargs='+', help='Target taxa in format taxon:value (e.g., order:coleoptera family:primates)')
    parser.add_argument('--log', help='Path to log file (if not provided, logs to stderr)', default=None)
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
                target_taxa[level.lower()] = value.lower()
            except ValueError:
                logger.warning(f"Invalid target taxa format '{item}'. Expected format is 'taxon:value'")
    
    main(tax_file=args.tax_results, fasta_file=args.fasta, lineage_file=args.rankedlineage_db, target_taxa=target_taxa)