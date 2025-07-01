#!/usr/bin/env python3

import sys
from collections import defaultdict

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
            # Positions for tax IDs: 6, 12, 18, 24
            # Positions for coverage: 9, 15, 21, 27
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

def main(file_path):
    """Main function to run the parser."""

    print(file_path)
    scaffold_taxa_lengths = parse_taxonomy_file(file_path)
    summary = generate_summary(scaffold_taxa_lengths)
    
    # Output the results
    print("scaffold_id,tax_id,length,top_n")
    for row in summary:
        print(f"{row[0]},{row[1]},{row[2]},{row[3]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python taxonomy_parser.py <taxonomy_file>")
        sys.exit(1)
    
    main(sys.argv[1])
