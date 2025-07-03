# FCS Taxonomy Parser

A simple Python script for parsing and summarizing taxonomy classification files from FCS (Foreign Contamination Screen) output.

## Usage

```bash
python3 fcs_taxonomy_parser.py -t <taxonomy_file> [-f <fasta_file>] [-r <rankedlineage_db>] [--target_taxa <taxon:value> ...]
```

### Arguments

- `-t, --tax_results`: Path to the taxonomy file to parse (required)
- `-f, --fasta`: Path to the FASTA file to process (optional, for calculating coverage percentages)
- `-r, --rankedlineage_db`: Path to NCBI ranked lineage database file (optional, for taxonomic filtering). [FTP link](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).
- `--target_taxa`: Target taxa in format taxon:value, taxon in lower-case (e.g., order:Coleoptera or family:Carabidae). List of taxons: ['genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

## Output Format

The script outputs a CSV format to stdout with the following columns:

- `scaffold_id`: The ID of the scaffold without coordinates
- `tax_id`: The taxonomic ID
- `covered_length`: The sum of coverage lengths for this taxonomic ID
- `top_n`: The rank of this taxonomic ID for the scaffold (1 = top match)
- `percentage_covered`: The percentage of the scaffold covered by this taxon (when FASTA file is provided)
- `is_target`: Boolean indicating whether the taxon matches specified target taxa (when lineage database is provided)

## Examples

### Basic usage:

```bash
./fcs_taxonomy_parser.py -t taxonomy_results.rpt > taxonomy_summary.csv
```

### Complete example with all options:

```bash
./fcs_taxonomy_parser.py -t taxonomy_results.rpt -f scaffolds.fasta -r rankedlineage.dmp --target_taxa phylum:chordata > taxonomy_summary.csv
```

The output will look like:

```
scaffold_id,tax_id,length,top_n,percentage_covered,is_target
HAP1_SCAFFOLD_1,9606,150000,1,75.00,True
HAP1_SCAFFOLD_1,10090,50000,2,25.00,True
HAP1_SCAFFOLD_2,9606,200000,1,100.00,True
```

This indicates that for `HAP1_SCAFFOLD_1`, the top taxonomic match is tax ID 9606 with a total coverage of 150000 (75% of the scaffold), followed by tax ID 10090 with 50000 coverage (25% of the scaffold). Both are marked as target taxa because they belong to phylum Chordata.
