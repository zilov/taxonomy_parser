# FCS Taxonomy Parser

A simple Python script for parsing and summarizing taxonomy classification files from FCS (Foreign Contamination Screen) output.

## Usage

```bash
python3 fcs_taxonomy_parser.py <taxonomy_file>
```

## Output Format

The script outputs a CSV format to stdout with the following columns:

- `scaffold_id`: The ID of the scaffold without coordinates
- `tax_id`: The taxonomic ID
- `length`: The sum of coverage lengths for this taxonomic ID
- `top_n`: The rank of this taxonomic ID for the scaffold (1 = top match)

## Example

```bash
./fcs_taxonomy_parser.py taxonomy_results.rpt > taxonomy_summary.csv
```

The output will look like:

```
scaffold_id,tax_id,length,top_n
HAP1_SCAFFOLD_1,9606,150000,1
HAP1_SCAFFOLD_1,10090,50000,2
HAP1_SCAFFOLD_2,9606,200000,1
```

This indicates that for `HAP1_SCAFFOLD_1`, the top taxonomic match is tax ID 9606 with a total coverage of 150000, followed by tax ID 10090 with 50000 coverage.
