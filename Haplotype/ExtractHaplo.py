import argparse
import os
import pandas as pd
import numpy as np
from io import StringIO

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Compute minor allele frequency (maf) from a VCF file.')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file path')
    parser.add_argument('-io', '--iodir', required=True, help='Input and Output directory for the CSV files')
    parser.add_argument('-fq', '--qtlfilename', required=True, help='Name for the QTL CSV file')
    parser.add_argument('-fp', '--panelfilename', required=True, help='Name for the panel CSV file')
    parser.add_argument('-fm', '--mergedfilename', required=True, help='Name for the merged CSV file')
    parser.add_argument('-q', '--qtl', required=True, help='QTL mapping file')
    parser.add_argument('-p', '--panel', required=True, help='Panel mapping file')
    parser.add_argument('-m', '--merged', required=True, help='Merged mapping file')

    return parser.parse_args()

def ensure_output_directory(output_dir):
    """Ensure the output directory exists."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def read_vcf(vcf_file):
    """Read the VCF file and return the data as a pandas DataFrame."""
    with open(vcf_file, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    
    if not lines:
        raise ValueError("The VCF file contains no data after headers.")
    
    return pd.read_csv(StringIO(''.join(lines)), sep='\t')

def validate_vcf_columns(vcf_df):
    """Ensure the VCF DataFrame contains the required columns."""
    required_columns = ['#CHROM', 'POS']
    if not all(col in vcf_df.columns[:9] for col in required_columns):
        raise ValueError("The VCF file is missing required columns.")

def separate_haplotypes(genotype_col):
    """Separate haplotypes from genotype strings."""
    genotype_col = genotype_col.astype(str)
    haplotypes = genotype_col.str.split('[|/]', expand=True, regex=True)
    # Convert haplotypes to numeric, coerce errors (e.g., invalid strings become NaN)
    haplotype1 = pd.to_numeric(haplotypes[0])
    haplotype2 = pd.to_numeric(haplotypes[1])

    return np.array([haplotype1.to_numpy(), haplotype2.to_numpy()])

def process_haplotypes(filtered_vcf_df):
    """Process haplotypes for all sample columns in the filtered VCF DataFrame."""
    haplo_list = []
    
    # Process each column starting from the 10th (sample data in VCF typically starts after metadata)
    for col_name in filtered_vcf_df.columns[9:]:
        haplo = separate_haplotypes(filtered_vcf_df[col_name])
        # Append the two haplotypes for the current sample column
        haplo_list.append(haplo[0])  # First haplotype
        haplo_list.append(haplo[1])  # Second haplotype
    
    # Combine the haplotypes into a 2D array for further analysis
    return np.array(haplo_list)

def main():
    args = parse_arguments()

    # Input and output paths
    vcf_file = args.input
    io_dir = args.iodir
    qtl_file = args.qtl
    panel_file = args.panel
    merged_file = args.merged
    Csv_name_snp = args.panelfilename
    Csv_name_qtl = args.qtlfilename
    Csv_name_mer = args.mergedfilename

    # Load mapping data
    input_file_qtl= os.path.join(io_dir, qtl_file)
    input_file_panel= os.path.join(io_dir, panel_file)
    input_file_mer= os.path.join(io_dir, merged_file)
    map_df_qtl = pd.read_csv(input_file_qtl, sep=",")
    map_df_panel = pd.read_csv(input_file_panel, sep=",")
    merged_df = pd.read_csv(input_file_mer, sep=",")

    ensure_output_directory(io_dir)

    # Read VCF data and validate
    vcf_df = read_vcf(vcf_file)
    validate_vcf_columns(vcf_df)

    merge_keys_left = ['#CHROM', 'POS']
    merge_keys_right = ['chr', 'pos']
    # Filter VCF data for SNP and QTL
    filtered_vcf_df_snp = pd.merge(vcf_df, map_df_panel[merge_keys_right], how='inner', left_on=merge_keys_left, right_on=merge_keys_right).drop(columns=merge_keys_right)
    filtered_vcf_df_qtl = pd.merge(vcf_df, map_df_qtl[merge_keys_right], how='inner', left_on=merge_keys_left, right_on=merge_keys_right).drop(columns=merge_keys_right)
    filtered_vcf_df_mer = pd.merge(vcf_df, merged_df[merge_keys_right], how='inner', left_on=merge_keys_left, right_on=merge_keys_right).drop(columns=merge_keys_right)
    # 排序 DataFrame
    filtered_vcf_df_snp = filtered_vcf_df_snp.sort_values(by=['#CHROM', 'POS'])
    filtered_vcf_df_qtl = filtered_vcf_df_qtl.sort_values(by=['#CHROM', 'POS'])
    filtered_vcf_df_mer = filtered_vcf_df_mer.sort_values(by=['#CHROM', 'POS'])

    # Process haplotypes for SNP and QTL
    combined_matrix_snp = process_haplotypes(filtered_vcf_df_snp)
    combined_matrix_qtl = process_haplotypes(filtered_vcf_df_qtl)
    combined_matrix_mer = process_haplotypes(filtered_vcf_df_mer)
    # Output CSV files
    output_file_snp = os.path.join(io_dir, Csv_name_snp)
    output_file_qtl = os.path.join(io_dir, Csv_name_qtl)
    output_file_mer = os.path.join(io_dir, Csv_name_mer)
    pd.DataFrame(combined_matrix_snp).to_csv(output_file_snp, index=False)
    pd.DataFrame(combined_matrix_qtl).to_csv(output_file_qtl, index=False)
    pd.DataFrame(combined_matrix_mer).to_csv(output_file_mer, index=False)

if __name__ == "__main__":
    main()

