#!/usr/bin/env python3

import sys
import argparse
import logging
from collections import defaultdict
from typing import Dict
import pandas as pd
import os

def setup_logger():
    """Configure and return a logger for the application."""
    logger = logging.getLogger('sourmash_taxonomy_parser')
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

def parse_sourmash_results(file_path):
    """Parse the sourmash results file and extract query-match relationships."""
    query_matches = defaultdict(list)
    
    try:
        df = pd.read_csv(file_path)
        logger.info(f"Loaded {len(df)} sourmash results")
        
        for _, row in df.iterrows():
            query_name = row['query_name']
            match_name = row['match_name']
            containment = row['containment']
            jaccard = row['jaccard']
            intersect_hashes = row.get('intersect_hashes', 0)  # Default to 0 if column doesn't exist

            if " " in match_name:
                # Handle cases where match_name contains spaces
                match_name = match_name.strip().split(" ")[0]
            
            query_matches[query_name].append({
                'match_name': match_name,
                'containment': containment,
                'jaccard': jaccard,
                'intersect_hashes': intersect_hashes
            })
        
        return query_matches
    except Exception as e:
        logger.error(f"Error parsing sourmash results file: {e}")
        raise

def parse_assembly_database(file_path):
    """Parse the assembly database file and return a DataFrame with taxonomic information."""
    try:
        assembly_df = pd.read_csv(file_path)
        logger.info(f"Loaded {len(assembly_df)} assembly records from database")
        
        # Set assembly name as index for faster lookups
        if 'assembly_accession' in assembly_df.columns:
            assembly_df.set_index('assembly_accession', inplace=True)
        else:
            # Use first column as index if standard names not found
            assembly_df.set_index(assembly_df.columns[0], inplace=True)
        
        return assembly_df
    except Exception as e:
        logger.error(f"Error parsing assembly database file: {e}")
        raise

def generate_summary(query_matches):
    """Generate summary of taxonomic matches for each query."""
    summary = []
    
    for query_name, matches in query_matches.items():
        # Sort matches by intersect_hashes in descending order (highest first)
        sorted_matches = sorted(matches, key=lambda x: x['intersect_hashes'], reverse=True)
        
        # Generate rows with rank
        for rank, match in enumerate(sorted_matches, 1):
            summary.append((
                query_name,
                match['match_name'],
                rank,
                match['containment'],
                match['jaccard'],
                match['intersect_hashes']
            ))
    
    return summary

def get_target_genomes(assembly_df, target_taxa):
    """Pre-compute set of all genomes that match target taxa."""
    if assembly_df is None or not target_taxa:
        return set()
    
    target_genomes = set()
    
    for level, target_value in target_taxa.items():
        if level in assembly_df.columns:
            # Vectorized search for this taxonomic level
            mask = (
                assembly_df[level].notna() & 
                (assembly_df[level].astype(str).str.lower().str.strip() == target_value) &
                (assembly_df[level].astype(str).str.strip() != 'nan')
            )
            level_matches = set(assembly_df[mask].index)
            target_genomes.update(level_matches)
            
            logger.info(f"Found {len(level_matches)} genomes matching {level}={target_value}")
    
    logger.info(f"Total target genomes: {len(target_genomes)}")
    return target_genomes

def write_summary_output(summary, output_file, assembly_df=None, target_taxa=None):
    """Write summary output to a file."""
    # Pre-compute target genomes set
    target_genomes = get_target_genomes(assembly_df, target_taxa)
    
    # Pre-compute taxid dictionary for faster lookup
    taxid_dict = {}
    if assembly_df is not None and 'taxid' in assembly_df.columns:
        taxid_dict = assembly_df['taxid'].dropna().astype(str).str.strip().to_dict()
    
    # Sort summary to ensure proper ordering
    sorted_summary = sorted(summary, key=lambda x: (x[0], -x[5]))
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("header,assembly_accession,taxa,top_n,containment,jaccard,intersect_hashes,is_target\n")
        
        for row in sorted_summary:
            query_name, match_name, rank, containment, jaccard, intersect_hashes = row
            
            # Fast lookup: is this genome in target set?
            is_target = match_name in target_genomes if target_genomes else None
            
            # Fast lookup: get taxid from pre-computed dictionary
            taxa = taxid_dict.get(match_name, "None")
            
            # Format the values
            is_target_str = str(is_target) if is_target is not None else "None"
            
            # Write the output
            f.write(f"{query_name},{match_name},{taxa},{rank},{containment:.6f},{jaccard:.6f},{intersect_hashes},{is_target_str}\n")

def write_non_target_output(summary_file, output_file, assembly_df):
    """Write non-target queries with lineage information of top1 match."""
    # Read the summary file to identify non-target queries
    df = pd.read_csv(summary_file)
    
    # Group by header and check if any row has is_target=True
    query_groups = df.groupby('header')
    non_target_queries = []
    
    for query_name, group in query_groups:
        # Check if this query has any target taxa matches
        has_target = group['is_target'].any()
        
        if not has_target:
            # Get the top1 match (top_n == 1)
            top1_rows = group[group['top_n'] == 1]
            if not top1_rows.empty:
                top1_row = top1_rows.iloc[0]
                non_target_queries.append(top1_row)
    
    # Write the non-target queries file
    with open(output_file, 'w') as f:
        # Write header
        f.write("header,assembly_accession,taxa,containment,jaccard,intersect_hashes,species,genus,family,order,class,phylum,kingdom\n")
        
        for _, row in enumerate(non_target_queries):
            query_name = row['header']
            match_name = row['assembly_accession']
            taxa = row['taxa']
            containment = row['containment']
            jaccard = row['jaccard']
            intersect_hashes = row['intersect_hashes']
            
            # Get lineage information
            lineage_info = ["None"] * 7  # species, genus, family, order, class, phylum, kingdom
            try:
                if match_name in assembly_df.index:
                    lineage_row = assembly_df.loc[match_name]
                    lineage_columns = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
                    for i, col in enumerate(lineage_columns):
                        if col in lineage_row.index and pd.notna(lineage_row[col]) and str(lineage_row[col]).strip() != 'nan':
                            lineage_info[i] = str(lineage_row[col]).strip()
            except Exception as e:
                logger.debug(f"Error getting lineage for {match_name}: {e}")
            
            # Write the output
            lineage_str = ",".join(lineage_info)
            f.write(f"{query_name},{match_name},{taxa},{containment:.6f},{jaccard:.6f},{intersect_hashes},{lineage_str}\n")

def main(sourmash_file=None, assembly_db=None, target_taxa=None, outdir=None):
    """Main function to run the sourmash parser."""
    
    # Initialize variables
    assembly_df = None
    query_matches = defaultdict(list)
    
    # Process assembly database if provided
    if assembly_db:
        logger.info(f"Processing assembly database: {assembly_db}")
        assembly_df = parse_assembly_database(assembly_db)
    
    if target_taxa:
        logger.info(f"Target taxa filters: {target_taxa}")
    
    logger.info(f"Processing sourmash results file: {sourmash_file}")
    query_matches = parse_sourmash_results(sourmash_file)
    
    # Generate summary
    summary = generate_summary(query_matches)
    
    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    
    # Generate output filenames
    sourmash_basename = os.path.basename(sourmash_file)
    summary_filename = f"{sourmash_basename}.summary.csv"
    summary_file = os.path.join(outdir, summary_filename)
    
    logger.info(f"Writing summary to: {summary_file}")
    
    # Write summary output
    write_summary_output(summary, summary_file, assembly_df, target_taxa)
    
    # Write non-target queries if assembly database and target taxa are available
    if assembly_df is not None and target_taxa:
        non_target_filename = f"{sourmash_basename}.non_target.csv"
        non_target_file = os.path.join(outdir, non_target_filename)
        logger.info(f"Writing non-target queries to: {non_target_file}")
        write_non_target_output(summary_file, non_target_file, assembly_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse sourmash results and extract query-match relationships with similarity scores.")
    parser.add_argument('-s', '--sourmash_results', help='Path to the sourmash results file to parse', required=True)
    parser.add_argument('-a', '--assembly_db', help='Path to assembly database file with taxonomic information', required=True)
    parser.add_argument('--target_taxa', nargs='+', help='Target taxa in format taxon:value (e.g., order:coleoptera family:Carabidae)')
    parser.add_argument('--log', help='Path to log file (if not provided, logs to stderr)', default=None)
    parser.add_argument('-o', '--outdir', help='Output directory for results (default: current directory)', default=None, required=True)
    args = parser.parse_args()
    
    # Configure file logging if requested
    if args.log:
        file_handler = logging.FileHandler(args.log)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
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
    
    main(sourmash_file=args.sourmash_results, assembly_db=args.assembly_db, target_taxa=target_taxa, outdir=args.outdir)
