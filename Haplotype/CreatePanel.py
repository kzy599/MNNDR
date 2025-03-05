import os
import argparse
import pandas as pd

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Compute minor allele frequency (maf) from a ChPoMa file.')
    parser.add_argument('-c', '--cpm', required=True, help='ChPoMa file')
    parser.add_argument('-io', '--iodir', required=True, help='Input and Output directory for the CSV files')
    parser.add_argument('-q', '--qtl', required=True, help='QTL mapping file')
    parser.add_argument('-p', '--panel', required=True, help='Panel mapping file')
    parser.add_argument('-m', '--merged', required=True, help='Merged mapping file')
    parser.add_argument('-s', '--panelsize', required=True, help='Panelsize')
    return parser.parse_args()

import pandas as pd
import numpy as np

def adjust_num_based_on_site(df_out, map_df, panel_size):
    """
    Adjusts the 'num' column in df_out based on the number of SNPs per chromosome in map_df.
    
    Args:
    - df_out: DataFrame with columns ['chr', 'len', 'num']
    - map_df: DataFrame with columns ['chr', 'pos', 'maf', 'QTL'] used to calculate the number of SNPs per chromosome
    - panel_size: The target total panel size for SNPs.
    
    Returns:
    - df_out: Adjusted DataFrame without the 'mod' column.
    """
    
    # Create numSite DataFrame similar to numSite in R
    num_site = map_df.groupby('chr').size().reset_index(name='num')
    
    # Initialize the mod column in df_out for tracking updates
    df_out['mod'] = True
    
    # First adjustment: update 'num' based on numSite
    for i, row in df_out.iterrows():
        matching_row = num_site[num_site['chr'] == row['chr']]
        if not matching_row.empty:
            matching_num = matching_row.iloc[0]['num']
            if matching_num < row['num']:
                df_out.at[i, 'num'] = matching_num
                df_out.at[i, 'mod'] = False

    # Loop while not all mod values are False (continue adjusting)
    while not all(df_out['mod'] == True):
        # Calculate the difference between target panel size and the current sum of 'num'
        cknumber = int(panel_size) - df_out['num'].sum()

        # Update the rows where 'mod' is True
        df_out.loc[df_out['mod'] == True, 'num'] += (
            df_out.loc[df_out['mod'] == True, 'len'] * cknumber / df_out.loc[df_out['mod'] == True, 'len'].sum()
        )

        # Round the updated 'num' values
        df_out.loc[df_out['mod'] == True, 'num'] = df_out.loc[df_out['mod'] == True, 'num'].round()

        # Adjust the largest 'num' to ensure the total sums match panel_size
        max_idx = df_out[df_out['mod'] == True]['num'].idxmax()
        df_out.at[max_idx, 'num'] += int(panel_size) - df_out['num'].sum()

        df_out['mod'] = True
        # Reapply the logic to update 'num' based on numSite again
        for i, row in df_out.iterrows():
            matching_row = num_site[num_site['chr'] == row['chr']]
            if not matching_row.empty:
                matching_num = matching_row.iloc[0]['num']
                if matching_num < row['num']:
                    df_out.at[i, 'num'] = matching_num
                    df_out.at[i, 'mod'] = False

    # Drop the 'mod' column as it's no longer needed
    df_out = df_out.drop(columns=['mod'])

    return df_out

def score_greedy(start, end, maf, pos):
    """Calculate the greedy score based on position and maf."""
    return maf * (end - start - abs((2 * pos) - end - start))

def greedy_choose_loci(num, map_df):
    """Greedy algorithm to select loci based on scores."""
    # Rename columns to match expected format
    num.columns = ['chr', 'length', 'num']
    map_df.columns = ['chr', 'pos', 'maf', 'QTL']
    panel = pd.DataFrame()

    # Iterate over each chromosome
    for i in range(len(num)):
        chr_name = num.loc[i, 'chr']
        cands = map_df[map_df['chr'] == chr_name].copy()

        if len(cands) < num.loc[i, 'num']:
            raise ValueError(f"Not enough SNPs in chromosome {chr_name}")
        
        # Initial score calculation for candidates
        cands['score'] = score_greedy(0, num.loc[i, 'length'], cands['maf'], cands['pos'])
        cands['lastStart'] = 0
        cands['lastEnd'] = num.loc[i, 'length']
        num_snps = 0
        
        while num_snps < num.loc[i, 'num']:
            if cands.empty:
                break

            # Select SNP with maximum score
            temp = cands['score'].idxmax()
            chosen = cands.loc[[temp]]
            panel = pd.concat([panel, chosen], ignore_index=True)
            cands = cands.drop(temp)
            num_snps += 1

            chosen_pos = chosen['pos'].values[0]

            # Update scores for remaining candidates
            to_update = cands[(cands['lastStart'] < chosen_pos) & (cands['lastEnd'] > chosen_pos)]
            intervals = to_update[['lastStart', 'lastEnd']].drop_duplicates()

            if intervals.empty:
                continue

            interval = intervals.iloc[0]
            last_start, last_end = interval['lastStart'], interval['lastEnd']

            # Update SNPs before and after the chosen position
            temp_before = (cands['lastStart'] == last_start) & (cands['pos'] < chosen_pos)
            cands.loc[temp_before, 'score'] = score_greedy(last_start, chosen_pos, cands.loc[temp_before, 'maf'], cands.loc[temp_before, 'pos'])
            cands.loc[temp_before, 'lastEnd'] = chosen_pos

            temp_after = (cands['lastStart'] == last_start) & (cands['pos'] > chosen_pos)
            cands.loc[temp_after, 'score'] = score_greedy(chosen_pos, last_end, cands.loc[temp_after, 'maf'], cands.loc[temp_after, 'pos'])
            cands.loc[temp_after, 'lastStart'] = chosen_pos

    return panel.sort_values(by=['chr', 'pos'])

def main():
    args = parse_arguments()

    # Input and output paths
    cpm_file = args.cpm
    io_dir = args.iodir
    qtl_file = args.qtl
    panel_file = args.panel
    merged_file = args.merged
    panel_size = args.panelsize
    # Read the map from ChPoMa file
    input_file_cpm = os.path.join(io_dir, cpm_file)
    map_df = pd.read_csv(input_file_cpm, sep=",")
    map_df.columns = ['chr', 'pos', 'maf', 'QTL']

    # Separate QTL and non-QTL SNPs
    map_df_qtl = map_df[map_df["QTL"] == True].copy().sort_values(by=['chr', 'pos'])
    map_df_snp = map_df[map_df["QTL"] != True].copy()

    # Define chromosome names and lengths
    num = pd.DataFrame({
        'chr': [
            "NC_088853.1", "NC_088854.1", "NC_088855.1", "NC_088856.1", "NC_088857.1",
            "NC_088858.1", "NC_088859.1", "NC_088860.1", "NC_088861.1", "NC_088862.1"
        ],
        'len': [
            76070991, 61469542, 61039741, 57946171, 57274926, 
            56905015, 53672946, 51133819, 50364239, 37310742
        ]
    })

    # Calculate the number of SNPs per chromosome based on lengths
    df_out = num.copy()
    df_out['num'] = (df_out['len'] * int(panel_size)) / df_out['len'].sum()
    df_out['num'] = df_out['num'].round()

    # Adjust SNP count to match total panel size
    difference = int(panel_size) - df_out['num'].sum()
    df_out.loc[df_out['num'].idxmax(), 'num'] += difference

    df_out = adjust_num_based_on_site(df_out=df_out, map_df=map_df_snp,panel_size=int(panel_size))
    # Run the greedy algorithm to select SNPs
    panel = greedy_choose_loci(num=df_out, map_df=map_df_snp)



    # Combine the two DataFrames row-wise
    panel_subset = panel[['chr', 'pos', 'maf', 'QTL']]
    merged_df = pd.concat([map_df_qtl, panel_subset], axis=0, ignore_index=True)
    merged_df.sort_values(by=['chr', 'pos'], inplace=True)

    # Write output files
    output_file_snp = os.path.join(io_dir, panel_file)
    output_file_qtl = os.path.join(io_dir, qtl_file)
    output_file_mer = os.path.join(io_dir, merged_file)

    panel.to_csv(output_file_snp, index=False)
    map_df_qtl.to_csv(output_file_qtl, index=False)
    merged_df.to_csv(output_file_mer, index=False)

if __name__ == "__main__":
    main()
