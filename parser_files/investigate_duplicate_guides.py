# This file works solely with the guide library file - does not change or alter the counts or read fastq files


# %% 0. Import/Setup

import pandas as pd
import numpy as np
import matplotlib.pyplot as plot
guides_file = "/Users/Claire/Downloads/git_clones/dual_guide_crispr_optimization/parser_files/20200513_library_1_2_unbalanced_dJR051.txt"

guides_df = pd.read_csv(guides_file, sep='\t')
# Check duplicate guides

# Convert all values in the 'Protospacer_A' column to uppercase
guides_df['protospacer_A'] = guides_df['protospacer_A'].str.upper()
# Convert all values in the 'Protospacer_B' column to uppercase
guides_df['protospacer_B'] = guides_df['protospacer_B'].str.upper()

# Subset and sort the cols of interest
guides_df_cols = guides_df[['gene', 'sgID_A', 'protospacer_A', 'sgID_B', 'protospacer_B']]

# %% Ignore
# Ignore lines below until next cell.
    # Lines originally written while learning python. Wrapped everything in lines into find_inconsistent function
# # Subset duplicate values in protospacer A col
#     # Subset and sort duplicated guide sequences in protospacer A column, keep=false specifies to mark all duplicates as TRUE
# all_dup_proto_a = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_A', keep=False)].sort_values(by='protospacer_A')
#     # Subset duplicated guide sequences in protospacer A column and return the unique values
# unique_of_dup_proto_a = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_A')]['protospacer_A'].unique()
#     # Subset the duplicates guides in protospacer_A column and return all instances that are truly duplicate values. First instance=FALSE, subsequent dups=TRUE,
# true_dup_proto_a = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_A')].sort_values(by='protospacer_A')
# print(f"There are {len(all_dup_proto_a)} total duplicate values within the protospacer A column, of these duplicate elements, {len(unique_of_dup_proto_a)} are unique and {len(true_dup_proto_a)} are true duplicates (first instance ignored)")
#     # Quick sanity check to ensure number printed is indeed the number of unique guides from duplicates
# proto_a_no_duplicates = all_dup_proto_a.drop_duplicates(subset='protospacer_A', keep='first')
# print(f"Sanity check that the number of {len(unique_of_dup_proto_a)} is equal to {len(proto_a_no_duplicates)}")
#
# # Subset duplicate values in protospacer B col
# all_dup_proto_b = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_B', keep=False)].sort_values(by='protospacer_B')
# unique_of_dup_proto_b = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_B')]['protospacer_B'].unique()
# true_dup_proto_b = guides_df_cols[guides_df_cols.duplicated(subset='protospacer_B')].sort_values(by='protospacer_B')
# print(f"There are {len(all_dup_proto_b)} total duplicate values within the protospacer B column, of these duplicate elements, {len(unique_of_dup_proto_b)} are unique and {len(true_dup_proto_b)} are true duplicates (first instance ignored)")
#     # Quick sanity check to ensure number printed is indeed the number of unique guides from duplicates
# proto_b_no_duplicates = all_dup_proto_b.drop_duplicates(subset='protospacer_B', keep='first')
# print(f"Sanity check that the number of {len(unique_of_dup_proto_b)} is equal to {len(proto_b_no_duplicates)}")
#
# # Ensure dups have the same gene
# # Group by the duplicate column 'A' and check the unique values in column 'B'
# grouped_a = all_dup_proto_a.groupby('protospacer_A')['gene'].nunique()
# grouped_b = all_dup_proto_b.groupby('protospacer_B')['gene'].nunique()
#
# # Display the groups where there are more than one unique value in column 'B'
# inconsistent_groups_a = grouped_a[grouped_a > 1]
# print(f"Inconsistent groups where a 'Protospacer A' corresponds with different values in 'sgID_A': {len(inconsistent_groups_a)}")
# # print(inconsistent_groups_a)
# inconsistent_groups_b = grouped_b[grouped_b > 1]
# print(f"Inconsistent groups where a 'Protospacer B' corresponds with different values in 'sgID_B': {len(inconsistent_groups_b)}")
# # print(inconsistent_groups_b)



# %% Create Find duplicates function
def find_inconsistent_guides(df, guide_column, gene_column):
    """
    Find rows with duplicate guide sequences where the gene names don't match.

    Args:
    df (pd.DataFrame): The input DataFrame.
    guide_column (str): The name of the guide sequence column.
    gene_column (str): The name of the gene column.

    Returns:
    pd.DataFrame: Rows with inconsistent gene names for the same guide sequence.
    """
    # Find duplicate guide sequences
    duplicates = df[df.duplicated(guide_column, keep=False)]

    # Group by guide sequence and check unique values in the gene column
    grouped = duplicates.groupby(guide_column)[gene_column].nunique()

    # Identify guide sequences with more than one unique gene
    inconsistent_guides = grouped[grouped > 1].index

    # Subset rows with inconsistent guide sequences
    inconsistent_rows = duplicates[duplicates[guide_column].isin(inconsistent_guides)]

    return inconsistent_rows

# Test dups in protospacer_A col
if __name__ == "__main__":
    # Find inconsistent guide sequences
    inconsistent_rows_a = find_inconsistent_guides(guides_df_cols, 'protospacer_A','gene').sort_values(by='protospacer_A')
    # Display the result
    print(f"\n Rows with inconsistent gene names for the same Protospacer A guide sequence: \n {inconsistent_rows_a}")

# Test dups in protospacer_B col
if __name__ == "__main__":
    # Find inconsistent guide sequences
    inconsistent_rows_b = find_inconsistent_guides(guides_df_cols, 'protospacer_B','gene').sort_values(by='protospacer_B')
    # Display the result
    print(f"Rows with inconsistent gene names for the same Protospacer B guide sequence: \n {inconsistent_rows_b}")

# %% 1. Guide A duplicates
## Look at duplicate sequences in protospacer A
inconsistent_rows_a['dup_seq_A'] = ''
inconsistent_rows_a['dup_seq_B'] = ''
inconsistent_rows_a['unique_seq_sgID_A'] = ''
inconsistent_rows_a['unique_seq_sgID_B'] = ''

# Identify Duplicate protospacer_A - create a boolean mask dup_A identifying rows where protospacer_A values are duplicated.
dup_A = inconsistent_rows_a['protospacer_A'].duplicated(keep=False)

# Process rows with duplicate 'protospacer_A'
for idx, row in inconsistent_rows_a[dup_A].iterrows():
    inconsistent_rows_a.at[idx, 'dup_seq_A'] = 'yes'
    #A subset of rows with the same protospacer_A value is created.
    subset = inconsistent_rows_a[inconsistent_rows_a['protospacer_A'] == row['protospacer_A']]
    #If any value in the protospacer_B column within this subset is duplicated, then examine sgID_B column
    if subset['protospacer_B'].duplicated(keep=False).any():
        inconsistent_rows_a.at[idx, 'dup_seq_B'] = 'yes'
        # Extract characters before the underscore in sgIDB and check for uniqueness
        sgIDB_parts = subset['sgID_B'].apply(lambda x: x.split('_')[0])
        if sgIDB_parts.nunique() > 1:
            inconsistent_rows_a.loc[(inconsistent_rows_a['protospacer_A'] == row['protospacer_A']) & (
                        inconsistent_rows_a['dup_seq_B'] == 'yes'), 'unique_seq_sgID_B'] = 'unique seq with ID B'

# Display the DataFrame
print(inconsistent_rows_a)



# %% 2. Guide B duplicates
## Look at duplicate sequences in protospacer A
inconsistent_rows_b['dup_seq_A'] = ''
inconsistent_rows_b['dup_seq_B'] = ''
inconsistent_rows_b['unique_seq_sgID_A'] = ''
inconsistent_rows_b['unique_seq_sgID_B'] = ''

# Identify Duplicate protospacer_A - create a boolean mask dup_A identifying rows where protospacer_A values are duplicated.
dup_B = inconsistent_rows_b['protospacer_B'].duplicated(keep=False)

# Process rows with duplicate 'protospacer_A'
for idx, row in inconsistent_rows_b[dup_B].iterrows():
    inconsistent_rows_b.at[idx, 'dup_seq_B'] = 'yes'
    #A subset of rows with the same protospacer_A value is created.
    subset = inconsistent_rows_b[inconsistent_rows_b['protospacer_B'] == row['protospacer_B']]
    #If any value in the protospacer_B column within this subset is duplicated, then examine sgID_B column
    if subset['protospacer_B'].duplicated(keep=False).any():
        inconsistent_rows_b.at[idx, 'dup_seq_B'] = 'yes'
        # Extract characters before the underscore in sgIDB and check for uniqueness
        sgIDA_parts = subset['sgID_A'].apply(lambda x: x.split('_')[0])
        if sgIDA_parts.nunique() > 1:
            inconsistent_rows_b.loc[(inconsistent_rows_b['protospacer_B'] == row['protospacer_B']) & (
                        inconsistent_rows_b['dup_seq_A'] == 'yes'), 'unique_seq_sgID_A'] = 'unique seq with ID B'

# Display the DataFrame
print(inconsistent_rows_b)


# TODO: mk df with duplicate values in proto A and proto B --> put into
# TODO: incoporate counts function from mageck to parser file
# TODO: generate summary tables that show total duplicates in A,
    # Maybe make ven diagram? - circle for A circle for B, overlaped area is duplicates in A and B
        # non overlaped area is duplicates in just A OR just B

# for duplicate values in proto a col, that also have dup values in proto b col, add counts to all sequnces




# %% Challenge
# Try rewriting it to take in the library guide df
# Look at duplicate sequences in protospacer A
guide_a_mismatch = guides_df_cols.copy()

guide_a_mismatch['dup_seq_A'] = ''
guide_a_mismatch['dup_seq_B'] = ''
guide_a_mismatch['unique_seq_sgID_A'] = ''
guide_a_mismatch['unique_seq_sgID_B'] = ''

# Identify duplicate values in 'protospacer_A'
dup_A = guide_a_mismatch['protospacer_A'].duplicated(keep=False)

# Use the function to find inconsistent guides
guide_a_mismatch = find_inconsistent_guides(guide_a_mismatch, 'protospacer_A', 'gene')

# Process rows with duplicate 'protospacer_A'
for idx, row in guide_a_mismatch[dup_A].iterrows():
    guide_a_mismatch.at[idx, 'dup_seq_A'] = 'yes'
    subset = guide_a_mismatch[guide_a_mismatch['protospacer_A'] == row['protospacer_A']]

    if subset['protospacer_B'].duplicated(keep=False).any():
        guide_a_mismatch.at[idx, 'dup_seq_B'] = 'yes'
        sgIDB_parts = subset['sgID_B'].apply(lambda x: x.split('_')[0])
        if sgIDB_parts.nunique() > 1:
            guide_a_mismatch.loc[(guide_a_mismatch['protospacer_A'] == row['protospacer_A']) & (
                    guide_a_mismatch['dup_seq_B'] == 'yes'), 'unique_seq_sgID_B'] = 'unique seq with ID B'

    # Label inconsistent guides
    if row['protospacer_A'] in inconsistent_guides:
        guide_a_mismatch.at[idx, 'unique_seq_sgID_A'] = 'inconsistent gene names'
# Display the DataFrame
print(guide_a_mismatch)


# %% Misc Sandbox

# unique_in_dup_r1 = guides_df[guides_df.duplicated(subset='protospacer_A')]['protospacer_A'].unique()
# duplicated_r1_keys = guides_df[guides_df.duplicated(subset='protospacer_A')]
# print(f"There are {len(unique_in_dup_r1)} unique gRNAs in the Protospacer A column, out of these, {len(duplicated_r1_keys)} are duplicate values")
# print(np.sort(duplicated_r1_keys)[::-1])
#
##guides_df[["gene", "sgID_A", "protospacer_A"]].sort_values(by="protospacer_A")

#
# duplicated_r2_keys = guides_df[guides_df.duplicated(subset='protospacer_B', keep=False)]['protospacer_B'].unique()
# print(len(duplicated_r2_keys), "guides are duplicated within the Protospacer B column of the guides library")
# duplicated_pairs = guides_df[(guides_df['protospacer_A'].duplicated(keep=False)) & (guides_df['protospacer_B'].duplicated(keep=False))][['protospacer_A', 'protospacer_B']].values.tolist()
# print(len(duplicated_pairs), "guides are duplicate values across Protospacer A column and the Protospacer B column in the guides library")
#

# def subset_duplicate_rows(dataframe, column_name):
#     """
#     Finds and returns rows with duplicate values in a specific column.
#
#     Args:
#     dataframe (pd.DataFrame): The input DataFrame.
#     column_name (str): The column to check for duplicates.
#
#     Returns:
#     pd.DataFrame: Rows containing duplicate values in the specified column.
#     """
#     # Find duplicate values in the specified column
#     duplicates = dataframe[dataframe.duplicated(column_name, keep=False)]
#     return duplicates
#
# # Example usage
# if __name__ == "__main__":
#     # Specify the column name to find duplicates
#     column_name = 'protospacer_A'
#     # Find duplicates
#     duplicate_rows = subset_duplicate_rows(guides_df, column_name)
#     # Display the duplicate rows
#     print("Rows with duplicate values in column '{}':".format(column_name))
#     print(duplicate_rows)
#
# duplicate_rows[["gene", "sgID_A", "protospacer_A"]].sort_values(by="protospacer_A")