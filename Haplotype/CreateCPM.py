import pandas as pd
import random
import argparse
import os
from io import StringIO

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Compute minor allele frequency (maf) from a VCF file.')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file path')
    parser.add_argument('-o', '--output', required=True, help='Output directory for the CSV file')
    parser.add_argument('-f', '--filename', required=True, help='Name for the CSV file')
    parser.add_argument('-n', '--nQtl', required=True, help='The number of qtls')
    args = parser.parse_args()

    # Input and output paths
    vcf_file = args.input
    output_dir = args.output
    Csv_name = args.filename

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the VCF file, skipping header lines that start with '##'
    with open(vcf_file, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]

    # Check if the file has any lines after skipping headers
    if not lines:
        raise ValueError("The VCF file contains no data after headers.")

    # Create a pandas DataFrame from the VCF data
    vcf_df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    # Ensure the DataFrame has the necessary columns
    required_columns = ['#CHROM', 'POS']
    if not all(col in vcf_df.columns[:9] for col in required_columns):
        raise ValueError("The VCF file is missing required columns.")

    # Sample 1000 random variants (rows) from the DataFrame
    num_variants = len(vcf_df)
    sample_size = min(1000, num_variants)  # Adjust if less than 1000 variants
    random.seed(42)  # For reproducibility; remove or change seed as needed
    qtls = random.sample(range(num_variants), sample_size)

    # Define a function to compute maf for each variant
    def compute_maf(row,is_qtl):
        # Extract genotype columns (from column 9 onwards)
        gt_cols = row.iloc[9:]
        # Ensure genotype columns are strings
        gt_cols = gt_cols.astype(str)
        # Extract alleles from genotype strings (e.g., '0|1')
        allele1 = gt_cols.str[0]
        allele2 = gt_cols.str[2]
        # Convert allele strings to numeric values, handle missing data
        allele1 = pd.to_numeric(allele1, errors='coerce')
        allele2 = pd.to_numeric(allele2, errors='coerce')
        # Combine alleles and drop NaN values
        alleles = pd.concat([allele1, allele2]).dropna()
        if len(alleles) == 0:
            maf = None
        else:
            # Compute minor allele frequency (maf)
            maf = alleles.sum() / (2 * len(allele1))
        # Return a Series with Chromosome, Position, and maf
        return pd.Series({'chr': row['#CHROM'], 'pos': row['POS'], 'maf': maf,'QTL':is_qtl})

    # Apply the compute_maf function to each sampled variant
    ChPoMa = vcf_df.apply(lambda row: compute_maf(row, row.name in qtls), axis=1)

    # Write the ChPoMa DataFrame to a CSV file in the output directory
    output_file = os.path.join(output_dir, Csv_name)
    ChPoMa.to_csv(output_file, index=False)

    # Optionally, print the DataFrame
    print("Computed minor allele frequencies for sampled variants:")
    print(ChPoMa)

if __name__ == "__main__":
    main()
