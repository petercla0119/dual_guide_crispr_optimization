# -*- coding: utf-8 -*-
"""dual_guide_parser_v10ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1WXHtAI4zVKt6PZPIYEvLqe4IpDI8Afb1

# dual guide parser - December 2021
 - **Project:** iNDI.
 - **Author(s):** Dan Ramos, Lirong Peng, Faraz Faghri, Mike Nalls, and Nicholas Johnson

---
### Quick Description:
- **Problem:** We need a method that preprocesses guides for experiments, something that parses fastqs to only include those guides included in R1 and R1 from the concensus guide list ... we also need to catalogue failed reads. We also need this to take MiSeq and other data.
- **Solution:** The workflow below sums it up pretty well. Let's test out some code on small iNDI datasets provided by Dan R. We have added support from SeqIO

### Workflow:
0.   Set up notebook.
1.   Import data, this includes concensus guides and R1 + R2 fastqs.
2.   Filter out poor quality reads.
3.   Detect whether any guide or read truncation is necessary.
4.   Identify matching read groups across R1 and R2.
5.   Reduce the datasets to the read groups that match in R1 and R2.
6.   Split into 'hits.\'  and 'recombinants.\'.
'hits.\' denotes read group matches and the protospacers match.
'recombinants.\' denotes read group matches but one or more protospacers does not.
7.   Export 'hits.\' and 'recombinants.\' per fastq.

### Notes on data for testing:
- **20200513_library_1_2_unbalanced_dJR051.csv** = All elements of the dual sgRNA library. Sequence from protospacer_A and protospacer_B columns must be present in the same row to be considered a match.
- **UDP0011_S5_R1_001.fastq.gz** = Final 19 bases of each read should match the final 19 bases of "protospacer_A" sequence from "20200513_library_1_2_unbalanced_dJR051.csv".
- **UDP0011_S5_R2_001.fastq.gz** = First 20 bases of each read should match reverse complement of "protospacer_B" sequence from "20200513_library_1_2_unbalanced_dJR051.csv".
** The major change in V3 of this code is matching sequences for the guides using all UPPER CASE bases intead of being case-sensitive."
** V8 automatically detects whether the first base of the reads and/or guides is G and acts accordingly. The result may be essentially running 19bp matches. Users may also set manually.
** V9 adds the ability to match guide 2 to read 1 and the reverse complement of guide 1 to read 2 in addition to the standard workflow.
** The above is also true for read 2.

# 0.   Set up notebook.
"""

# %% Set up notebook

import os

import os
# from google.colab import drive
import numpy as np
import pandas as pd
import math
import sys
# import joblib
import subprocess
import argparse
import gzip
import textwrap
import warnings

#!pip install --upgrade tables
#! pip install biopython

# import tables
# import requests
# from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import islice
# Comment out below after testing.

#drive.mount('/content/drive/')
#os.chdir("/content/drive/Shared drives/CARD_iNDI/scratch/dual_guide_parser")
# ! pwd

# Set  options for testing.

# guides_file = "/Users/Claire/Downloads/git_clones/dual_guide_crispr_optimization/parser_files/20200513_library_1_2_unbalanced_dJR051.txt"
# r1_file = "/Users/Claire/Downloads/raw_sequencing/JH8105_1_S1_L001_R1_001.fastq.gz"
# r2_file = "/Users/Claire/Downloads/raw_sequencing/JH8105_1_S1_L001_R2_001.fastq.gz"
# N_rows = 1500 # Speed up testing, this just reads the first 10K sequences.
# check_length = 500 # top and bottom of the array, how far to check for whether composed with G.
# guide_1_offset = -999 # -999 is the sentinel value
# guide_2_offset = -999 # -999 is the sentinel value
# read_1_offset = -999 # -999 is the sentinel value
# read_2_offset = -999 # -999 is the sentinel value
# purity = 0.95
# check_reverse = True

# Set the options for production.

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

#Thanks for trying the dual_guide_parser from CARD + iNDI + DTi.
#To run this code, you will need to specify the guides_file, this is a file similar to the example 20200513_library_1_2_unbalanced_dJR051.csv.
#You will need to specify a pair of R1 and R2 files, such as UDP0007_S1_R1_001.fastq.gz and UDP0007_S1_R2_001.fastq.gz.
#You can also specify the number of read groups you are interested for testing the tool, this relates to the option n_groups.
#This will only read that many readgroups from the R1 and R2 files, allowing you to speed things up a bit.
#This code must be run from the working directory that contains the R1 and R2 files, but the guides_file can be anywhere, just
#specify a full path to the guides file like ~/Desktop/20200513_library_1_2_unbalanced_dJR051.csv.
#This code attempts to automatically equalize guides and reads. In other words,
#it will attempt to cut out letters that begin a guide as a rule, or a read. You may also specify these settings
#manually.
#This is best run on a large RAM / high CPU set up as the files are quite large.
#Finally, to run this code, you will need several packages, including biopython. To see the required packages listed, run with the -h option.
#
#'''))
parser.add_argument('--packages', help='Request for packages required to run, and how to install.', action='store_true')
parser.add_argument('--guides_file', type=str, default='missing', help='Mandatory input filepath. This is a file similar to the example 20200513_library_1_2_unbalanced_dJR051.csv. This can be a complete filepath')
parser.add_argument('--r1_file', type=str, default='missing', help='Mandatory input file name. An R1 file in your working directory.')
parser.add_argument('--r2_file', type=str, default='missing', help='Mandatory input file name. An R2 file in your working directory.')
parser.add_argument('--N_reads', type=int, default=0, help='Optional number of readgroups to test. An integer.')
parser.add_argument('--guide_1_offset', type=int, default = -999, help='# of characters to truncate from the left of the protospacer_A. Read truncation is not automatic, set accordingly. Default = 0.')
parser.add_argument('--read_1_offset', type=int, default = -999, help='# of characters to truncate from the left of read_1. Read truncation is not automatic, set accordingly. Default = 0.')
parser.add_argument('--guide_2_offset', type=int, default = -999, help='# of characters to truncate from the left of the protospacer_B. A Read truncation is not automatic, set accordingly. Default = 0.')
parser.add_argument('--read_2_offset', type=int, default = -999, help='# of characters to truncate from the left of read_2. Read truncation is not automatic, set accordingly. Default = 0.')
parser.add_argument('--purity', type=float, default = 0.95, help='minimum percentage of reads with the same beginning to be cut.')
parser.add_argument('--check_length', type=int, default=500, help='# of rows on the bottom and top of the array to check for concordance of starting letter.')
parser.add_argument('--check_reverse', action='store_true', help='# match reads against both guide1+reverse_comp(guide2) and guide2_reverse_comp(guide1)')

args = parser.parse_args()

if(args.packages):
  print("Must have numpy and pandas available. \n Additionally, must install biopython.")
  print("To install biopython, run pip install biopython, or, if using conda,")
  print("conda install -c conda-forge biopython")
  quit()

print("#"*46)
print("")
print("Here is some basic info on the command you are about to run.")
print("Python version info...")
print(sys.version)
print("CLI argument info...")
print("The guides file you are using is", args.guides_file, ".")
print("The r1 file you are using is", args.r1_file, ".")
print("The r2 file you are using is", args.r2_file, ".")
print("How many read groups are only for a quick test and not the full set?", args.N_reads, ".")
print("")
print("#"*46)

guides_file = args.guides_file
r1_file = args.r1_file
r2_file = args.r2_file
N_reads = args.N_reads
N_rows = args.N_reads
guide_1_offset = args.guide_1_offset
guide_2_offset = args.guide_2_offset
read_1_offset = args.read_1_offset
read_2_offset = args.read_2_offset
check_length = args.check_length
purity = args.purity
check_reverse = args.check_reverse

# Commented out IPython magic to ensure Python compatibility.
# %tb
# %% 1. Import data
"""# 1. Import data, this includes concensus guides and R1 + R2 fastqs."""

# Function to read the guide library file in .csv and .txt format
# TODO: Ensure function runs
# def read_guides_file(guides_file):
#     file_extension = os.path.splitext(guides_file)[1].lower()
#
#     if file_extension == '.txt':
#         # Read as a tab-delimited file
#         guides_df = pd.read_csv(guides_file, sep = '\t', engine = 'c')
#     elif file_extension == '.csv':
#         # Read as a comma-separated file
#         guides_df = pd.read_csv(guides_file, engine = 'c')
#     else:
#         raise ValueError("Unsupported file type: {}".format(file_extension))
#
#     return guides_df
#
# guides_df = read_guides_file(guides_file)

# Modify based on guide file format - Guide file received was in txt format
guides_df = pd.read_csv(guides_file, sep='\t')

"""# # Import R1s and R2s.
# ## pysam way Pysam introduces many other file format compatibilities such as
# ## CRAM/SAM/BAM
# if (N_rows == 0):
#   with pysam.FastxFile(r1_file) as fh:
#     r1_df = pd.DataFrame([(entry.name, entry.sequence, entry.comment, \
#     entry.quality) for entry in fh], columns=['name','seq', 'comment', 'qual'])
#   with pysam.FastxFile(r2_file) as fh:
#     r2_df = pd.DataFrame([(entry.name, entry.sequence, entry.comment, \
#     entry.quality) for entry in fh], columns=['name','seq', 'comment', 'qual'])
# else:
#   with pysam.FastxFile(r1_file) as fh:
#     r1_df = pd.DataFrame([(entry.name, entry.sequence, entry.comment, \
#     entry.quality) for entry in islice(fh,0,N_rows)],
#            columns=['name','seq', 'comment', 'qual'])
#   with pysam.FastxFile(r2_file) as fh:
#     r2_df = pd.DataFrame([(entry.name, entry.sequence, entry.comment, \
#     entry.quality) for entry in islice(fh,0,N_rows)],
#     columns=['name','seq', 'comment', 'qual'])
#
# ## SeqIO way
# For more compatitibility with other files types will need to import as SeqIO objects.
# Consider the following suggestions if so.
# use to_dict for compatibility with more file types and for true dictionary
# functionality. If files too big, use .index.
# Otherwise, if more memory needed, instantiate as a list """

# With open the gzip files (r1 and r2) in reading mode ('rt') and assign it to variable ('r1' and 'r2', respectively).
# Initialize 'r1_it' iterator to iterate over contents in the first file with 'FastqqGeneralIterator'
# Use of backslash ('\') is a continuation of the 'with' statement across multiple lines of code for better readability
# Check 'N_rows' is empty. Then create new dataframe from iterator. Dataframe has three columns.
# Slice 'r1_it' and 'r2_it' from (0) to 'N_rows'. Then create dataframes from the sliced iterator
with gzip.open(r1_file, mode = 'rt') as r1, \
     gzip.open(r2_file, mode = 'rt') as r2:
  r1_it = FastqGeneralIterator(r1)
  r2_it = FastqGeneralIterator(r2)

  if (N_rows == 0):
    r1_df = pd.DataFrame(r1_it, columns=['title', 'seq', 'qual'])
    r2_df = pd.DataFrame(r2_it, columns=['title', 'seq', 'qual'])
  else:
    r1_df = pd.DataFrame(islice(r1_it, 0, N_rows),
                         columns=['title', 'seq', 'qual'])
    r2_df = pd.DataFrame(islice(r2_it, 0, N_rows),
                         columns=['title', 'seq', 'qual'])

# Get the index of all sequences with at least one unacceptable quality base
# Access 'seq' column in 'r1_df' and check each string to see it contains the character 'N'. Returns bool
removed_r1 = r1_df.seq.str.contains('N')
removed_r2 = r2_df.seq.str.contains('N')

all_removed = removed_r1 | removed_r2

r1_df = r1_df.loc[~all_removed,:]
r2_df = r2_df.loc[~all_removed,:]

# Check the composition of the guides and the reads. Do they start with G?
# How long are the reads?
guide_1_length = len(guides_df.iloc[0,guides_df.columns.get_loc('protospacer_A')])
guide_2_length = len(guides_df.iloc[0,guides_df.columns.get_loc('protospacer_B')])
original_read_1_length = len(r1_df.iloc[0,r1_df.columns.get_loc('seq')])
original_read_2_length = len(r2_df.iloc[0,r1_df.columns.get_loc('seq')])

assert original_read_1_length >= guide_1_length, 'We must assume read 1 is at least as large as the guide. Otherwise, we did something wrong'
assert original_read_2_length >= guide_2_length, 'We must assume read 2 is at least as large as the guide. Otherwise, we did something wrong'

# Possible output for later, or usage here. The first letter composition of each thing
guide_1_frst_ltrs = guides_df.loc[0:check_length,'protospacer_A'].astype(str).str[0].value_counts()
guide_2_frst_ltrs = guides_df.loc[0:check_length,'protospacer_B'].astype(str).str[0].value_counts()
r1_frst_ltrs = r1_df.loc[0:check_length,'seq'].astype(str).str[0].value_counts()
r2_frst_ltrs = r2_df.loc[0:check_length,'seq'].astype(str).str[0].value_counts()

print("Guide and read truncation readout: ##########")
print(f"Length of input guide 1: {guide_1_length}")
print(f"Length of input guide 2: {guide_2_length}")
print(f"Length of input read 1: {original_read_1_length}")
print(f"Length of input read 2: {original_read_2_length}")

guide_1_end = guide_1_length
guide_2_end = guide_2_length

def get_offset(df, c_length, purity, column, ltr):
  """
  Get the offset based on the inputs.

  Parameters
  df - the dataframe to examine
  c_length - the number of rows to check for the letter
  purity - minimum percentage of reads that start with ltr
  column - the column to look in
  ltr - the starting letter to check for

  Returns
  offset_value - will be either 0 or 1
  """
  first_gs = df.iloc[:c_length, df.columns.get_loc(column)]
  last_gs = df.iloc[-c_length:, df.columns.get_loc(column)]
  # Check whether all rows start with 'G'
  first = sum(first_gs.str.startswith(ltr))
  last = sum(last_gs.str.startswith(ltr))
  if (((first + last) / (2 * c_length)) >= purity):
    offset_value = 1
  else:
    offset_value = 0
  return offset_value

# If the user did not guide_1/read_1 offsets
if guide_1_offset == -999: # -999 is the sentinel value meaning unset
  # If the guide_1 all start with 'G'
  guide_1_offset = get_offset(guides_df, check_length, purity, 'protospacer_A', 'G')
  guide_1_length = (guide_1_end - guide_1_offset)

if guide_2_offset == -999: # -999 is the sentinel value meaning unset
  # If the guide_2 all start with 'G'
  guide_2_offset = get_offset(guides_df, check_length, purity, 'protospacer_B', 'G')
  guide_2_length = guide_2_end - guide_2_offset

if read_1_offset == -999: # -999 is the sentinel value meaning unset
  # If the read_1 all start with 'G'
  read_1_offset = get_offset(r1_df, check_length, purity, 'seq', 'G')

if read_2_offset == -999: # -999 is the sentinel value meaning unset
  # If the read_2 all start with 'G'
  read_2_offset = get_offset(r2_df, check_length, purity, 'seq', 'C')

# Check everything
if read_1_offset + guide_1_length > original_read_1_length:
  raise IndexError("Read 1 not long enough to sustain truncation. " + \
                   "Check input and parameters or manually set offsets.")
if read_2_offset + guide_2_length > original_read_2_length:
  raise IndexError("Read 2 not long enough to sustain truncation " + \
                   "Check input and parameters or manually set offsets.")

# set read lengths based on guide lengths. As above, if the reads are too short
# there is something wrong with the input and things should be rethought.
# As long as the read is longer than the guide, the below operations should be
# safe.

read_1_length = guide_1_length
read_1_end = read_1_offset + read_1_length

read_2_length = guide_2_length
read_2_end = read_2_offset + read_2_length

if read_1_end - read_1_offset > guide_1_length:
  raise IndexError(f"After truncation, read 1 length ({str(read_1_length)})" + \
                  f"is too big! guide 1 length ({str(guide_1_length)})!")
if read_2_end - read_2_offset > guide_2_length:
  raise IndexError(f"After truncation, read 2 length ({str(read_2_length)})" + \
                    f"is too big! guide 2 length ({str(guide_2_length)})!")

print(f"Therefore, guide 1 is now {guide_1_length} bp and guide 2 is {guide_2_length} bp")
print(f"guide_1: Going from {guide_1_offset} to {guide_1_end} and guide_2:" + \
      f" from {guide_2_offset} to {guide_2_end}")
print(f"read_1: Going from {read_1_offset} to {read_1_end} and read_2" + \
      f" from {read_2_offset} to {read_2_end}")
print("The reads have been truncated to be the same size as the guides.")

assert read_1_length == guide_1_length, 'Oh no, fix read_1 and g1!'
assert read_2_length == guide_2_length, 'Oh no, fix read_2 and g2!'
assert read_1_end - read_1_offset == read_1_length, 'Oh no, fix read_1 and g2!'
assert read_2_end - read_2_offset == read_2_length, 'Oh no, fix read_2 and g2!'
assert guide_1_length == guide_1_end - guide_1_offset, 'Oh no, fix guide_1!'
assert guide_2_length == guide_2_end - guide_2_offset, 'Oh no, fix guide_2!'

# only do this once
r1_df.insert(loc=2, column='plus', value='+')
r2_df.insert(loc=2, column='plus', value='+')

def split_str(s: str) -> str:
    return s.split(" ", maxsplit = 1)[0]
r1_df['read_group'] = r1_df["title"].apply(split_str)
r2_df['read_group'] = r2_df["title"].apply(split_str)

r1_df['title'] = '@' + r1_df['title']
r2_df['title'] = '@' + r2_df['title']

"""This next code block defines the r1 and r2 keys that are often 19 BP, also making a reverse complement of r2 key."""

# New as of v5.
guides_df['protospacer_A_19bp_trimmed'] = [x[guide_1_offset:guide_1_end] \
                                           for x in guides_df['protospacer_A']]
guides_df['protospacer_B_19bp_trimmed'] = [x[guide_2_offset:guide_2_end] \
                                           for x in guides_df['protospacer_B']]

# Make guide key columns.
guides_df['r1_key'] = guides_df['protospacer_A_19bp_trimmed']

# R2 is tricky as it is the reverse compliment. Flip it and translate using the function below.
guides_df['r2_key'] = guides_df['protospacer_B_19bp_trimmed'].apply(reverse_complement)
guides_df['r1_r2_key'] = guides_df['r1_key'] + "_" + guides_df['r2_key']

# Include an option to do the reverse as well
if check_reverse:
  guides_df['r2_r1_key'] = guides_df['protospacer_B_19bp_trimmed'] + "_" + guides_df['protospacer_A_19bp_trimmed'].apply(reverse_complement)

# Get the guide seqs, these relate to the 19 BP segments in the guide_df.
r1_df.loc[:,'guide_seq'] = [x[read_1_offset:read_1_end] for x in r1_df.seq]
r2_df.loc[:,'guide_seq'] = [x[read_2_offset:read_2_end] for x in r2_df.seq]

# %% 2. Identify matching read groups

"""# 2. Identify matching read groups across R1 and R2."""

# Pull the read groups from R1 and R2. Make a concensus read group list.
r1_read_groups_df = r1_df[['read_group']]
r2_read_groups_df = r2_df[['read_group']]

consensus_read_groups_df = r1_read_groups_df.merge(r2_read_groups_df, on ='read_group', how='inner')

# Quantify obvious pair failures.

r1_N_attempted_read_groups = r1_df.read_group.shape[0]
r2_N_attempted_read_groups = r1_df.read_group.shape[0]
consensus_N_read_groups = consensus_read_groups_df.shape[0]
r1_pcnt_consensus = round((consensus_N_read_groups/r1_N_attempted_read_groups)*100, 2)
r2_pcnt_consensus = round((consensus_N_read_groups/r2_N_attempted_read_groups)*100, 2)

print("#"*46)
print(f"R1 had {r1_N_attempted_read_groups} potential read groups, of these {r1_pcnt_consensus} % were among the concensus read groups. \n"
      f"R2 had {r2_N_attempted_read_groups} potential read groups, of these {r2_pcnt_consensus} % were among the concensus read groups. \n"
      f"In total there are {consensus_N_read_groups} read groups that are matching across R1 and R2 for this experiment.")
print(f"{sum(all_removed)} potential read groups were removed for poor quality reads")
print("#"*46)

# %% 3. Reduce R1 and R2 to only those in the concensus_read_groups df.
"""# 3. Reduce the datasets to the read groups that match in R1 and R2.
"""

# Reduce R1 and R2 to only those in the concensus_read_groups df.
consensus_read_list = consensus_read_groups_df.read_group.unique()

r1_df['in_consensus'] = r1_df.read_group.isin(consensus_read_list)
r2_df['in_consensus'] = r2_df.read_group.isin(consensus_read_list)

r1_reduced_df = r1_df[r1_df['in_consensus'] == True]
r2_reduced_df = r2_df[r1_df['in_consensus'] == True]

# %% 4. Split data into two dataframes
"""# 4. Split into 'hits.\'  and 'recombinants.\'.  
'hits.\' denotes read group matches and the protospacers match.
'recombinants.\' denotes read group matches but one or more protospacers does not.
"""

# Build dataset that is used to check for recombination w/in each read group.
r1_guide_read_df = r1_reduced_df[['read_group', 'guide_seq']].copy()
r1_guide_read_df.rename(columns={'guide_seq': 'r1_guide_seq'}, inplace=True)

r2_guide_read_df = r2_reduced_df[['read_group', 'guide_seq']].copy()
r2_guide_read_df.rename(columns={'guide_seq': 'r2_guide_seq'}, inplace=True)

# # Build dataset that is used to check for recombination w/in each read group.
# r1_guide_read_df = r1_reduced_df[['read_group', 'guide_seq']]
# r1_guide_read_df.rename(columns={'guide_seq':'r1_guide_seq'}, inplace=True)
#
# r2_guide_read_df = r2_reduced_df[['read_group', 'guide_seq']]
# r2_guide_read_df.rename(columns={'guide_seq':'r2_guide_seq'}, inplace=True)

combined_guide_read_df = r1_guide_read_df.merge(r2_guide_read_df, on='read_group', how='inner')
combined_guide_read_df['combined_guide_seqs'] = combined_guide_read_df['r1_guide_seq'] + "_" + combined_guide_read_df['r2_guide_seq']

# Flag expected pairs from the guides_df.
reference_list = guides_df['r1_r2_key'].tolist()

if check_reverse:
  reference_list.extend(guides_df['r2_r1_key'].tolist())

uppercase_reference_list = [x.upper() for x in reference_list]

combined_guide_read_df['uppercase_combined_guide_seqs'] = combined_guide_read_df['combined_guide_seqs'].str.upper()

combined_guide_read_df['non_recombinant'] = combined_guide_read_df['uppercase_combined_guide_seqs'].isin(uppercase_reference_list)

# Counting hits
try:
    on_target = combined_guide_read_df['non_recombinant'].value_counts()[True]
except KeyError:
    warnings.warn("WARNING: There are no on-target hits. Something is probably wrong.",
                  category = RuntimeWarning)
    on_target = 0

# Counting recombinants
try:
    recombinant = combined_guide_read_df['non_recombinant'].value_counts()[False]
except KeyError:
    warnings.warn("WARNING: There are no recombinant hits. Something is probably wrong.",
                  category = RuntimeWarning)
    recombinant = 0

total_reads = on_target + recombinant
on_target_pcnt = round((on_target/total_reads)*100, 2)
recombinant_pcnt = round((recombinant/total_reads)*100, 2)

print("#"*46)
print(f"There are a total of {total_reads} potential read groups after filtering, of these {on_target_pcnt} % were on target for R1 and R2. This means {recombinant_pcnt} % are recombinant read groups.")
print("#"*46)
print('\n')
if(on_target_pcnt < 25):
  warnings.warn(f"Not many hits. Check if guides and reads were trimmed " + \
                f"appropriately. \nread_1 comp: \n{r1_frst_ltrs} \n" + \
                f"read_2 comp: \n{r2_frst_ltrs} \n " + \
                f"guide_1 comp: \n{guide_1_frst_ltrs} \n " + \
                f"guide_2 comp: \n{guide_2_frst_ltrs} ")

# Split data into recombinant and hit subsets.
# Select read groups where combined_guide_read_df['non_recombinant'] == True
hits_df = combined_guide_read_df[combined_guide_read_df['non_recombinant'] == True][['read_group']]
# Select read groups where combined_guide_read_df['non_recombinant'] == False
recombinant_df = combined_guide_read_df[combined_guide_read_df['non_recombinant'] == False][['read_group']]

hits_list = hits_df.read_group.tolist()

r1_reduced_df['hit'] = r1_reduced_df['read_group'].isin(hits_list)
r2_reduced_df['hit'] = r2_reduced_df['read_group'].isin(hits_list)

r1_hits_df = r1_reduced_df[r1_reduced_df['hit'] == True]
r1_recombinant_df = r1_reduced_df[r1_reduced_df['hit'] == False]

r2_hits_df = r2_reduced_df[r2_reduced_df['hit'] == True]
r2_recombinant_df = r2_reduced_df[r2_reduced_df['hit'] == False]

# %% 5. Export two dfs per fastq
"""# 5. Export 'hits.\' and 'recombinants.\' per fastq.
"""

# Now its just back to the fastqs from here ... ouch. Start stacking the hits!
r1_hits_stacked_df = r1_hits_df[['title', 'seq', 'plus', 'qual']].stack()
r2_hits_stacked_df = r2_hits_df[['title', 'seq', 'plus', 'qual']].stack()

r1_hits_fastq_df = r1_hits_stacked_df
r2_hits_fastq_df = r2_hits_stacked_df

# Identify fails, these are recombinants with either R1 or R2 not in the guide list. If the uppercase guide sequences are not in a match to uppercase r1_key or r2_key from guides_df.
uppercase_r1_keys = [x.upper() for x in guides_df['r1_key']]
uppercase_r2_keys = [x.upper() for x in guides_df['r2_key']]

r1_recombinant_df = r1_recombinant_df.copy()  # Ensure we're working with a copy
r1_recombinant_df['uppercase_guide_seq'] = r1_recombinant_df['guide_seq'].str.upper()
r2_recombinant_df = r2_recombinant_df.copy()  # Ensure we're working with a copy
r2_recombinant_df['uppercase_guide_seq'] = r2_recombinant_df['guide_seq'].str.upper()
# r1_recombinant_df.loc[:, 'uppercase_guide_seq'] = [x.upper() for x in r1_recombinant_df.loc[:, 'guide_seq']]
# r2_recombinant_df.loc[:, 'uppercase_guide_seq'] = [x.upper() for x in r2_recombinant_df.loc[:, 'guide_seq']]

r1_recombinant_df = r1_recombinant_df.copy()  # Ensure it's a copy, not a view
r1_recombinant_df.loc[:, 'in_guide_library'] = r1_recombinant_df['uppercase_guide_seq'].isin(uppercase_r1_keys)
r2_recombinant_df = r1_recombinant_df.copy()  # Ensure it's a copy, not a view
r2_recombinant_df.loc[:, 'in_guide_library'] = r2_recombinant_df['uppercase_guide_seq'].isin(uppercase_r2_keys)
# r1_recombinant_df.loc[:, 'in_guide_library'] = r1_recombinant_df.loc[:, 'uppercase_guide_seq'].isin(uppercase_r1_keys)
# r2_recombinant_df.loc[:, 'in_guide_library'] = r2_recombinant_df.loc[:, 'uppercase_guide_seq'].isin(uppercase_r2_keys)

# TODO: consider adding function for user input to select permissiveness of output files
    # TODO: matches r1 sequence
    # TODO: matches r2 sequence

"""# Split R1 and R2 into true recombinants and failed recombinants based on the 'in_guide_library' flag. But first make a list of read groups that fail. Then pull these readgroups to split true recombinants and fails.
"""
# May be due to recombination or sequence read error, need to differentate between the two reasons
r1_failed_recombinants = r1_recombinant_df[r1_recombinant_df['in_guide_library'] == False]
r2_failed_recombinants = r2_recombinant_df[r2_recombinant_df['in_guide_library'] == False]

# TODO: add function to calculate hamming distance and tolerate nucleotide substitution of a user-specified value

recombinant_failed_readgroups = pd.concat((r1_failed_recombinants['read_group'],
                                          r2_failed_recombinants['read_group']),
                                          axis=0
                                          ).unique()

N_big_fails = len(recombinant_failed_readgroups)
# FIXME: The sum of reads not in the guide library(len(r1_failed_recombinants) + len(r2_failed_recombinants)  should equal to N_big_fails = len(recombinant_failed_readgroups) since those sequences were not in the guide library
    # remember that uppercase_r2_keys maps to the reverse complement of protospacer B
print(f'Number of R1 failed recombinants: {len(r1_failed_recombinants)}')
print(f'Number of R2 failed recombinants: {len(r2_failed_recombinants)}')
print(f'The sum of R1 and R2 failed recombinants: {len(r1_failed_recombinants) + len(r2_failed_recombinants)}')
print(f'Number of failed read groups: {N_big_fails}')

# Creates a new column to flag if read groups are in 'recombinant_failed_readgroups'. If yes, then 'big_fail' == True, otherwise False
r1_recombinant_df = r1_recombinant_df.copy()  # Ensure it's a copy, not a view
r1_recombinant_df.loc[:, 'big_fail'] = r1_recombinant_df.loc[:, 'read_group'].isin(recombinant_failed_readgroups)
r2_recombinant_df = r2_recombinant_df.copy()  # Ensure it's a copy, not a view
r2_recombinant_df.loc[:, 'big_fail'] = r2_recombinant_df.loc[:, 'read_group'].isin(recombinant_failed_readgroups)
# r1_recombinant_df.loc[:, 'big_fail'] = r1_recombinant_df.loc[:, 'read_group'].isin(recombinant_failed_readgroups)
# r2_recombinant_df.loc[:, 'big_fail'] = r2_recombinant_df.loc[:, 'read_group'].isin(recombinant_failed_readgroups)
# FIXME: weird... there are 58 rows in r1_recombinant_df that are in the guide library AND a big fail... how?

# Read groups/guide seq are in 'recombinant_failed_readgroups'... TODO: Consider sequencing base error?
r1_failed_recombinants_df = r1_recombinant_df[r1_recombinant_df.loc[:, 'big_fail'] == True]
r2_failed_recombinants_df = r2_recombinant_df[r2_recombinant_df.loc[:, 'big_fail'] == True]

# print(f'{len(r1_recombinant_df[r1_recombinant_df['in_guide_library']==True])}')
# print(f'{len(r2_recombinant_df[r2_recombinant_df['in_guide_library']==True])}')
# print('\nNot in guide library')
# print(f'{len(r1_recombinant_df[r1_recombinant_df['in_guide_library']==False])}')
# print(f'{len(r2_recombinant_df[r2_recombinant_df['in_guide_library']==False])}')

# Read groups/guide seq are
r1_true_recombinants_df = r1_recombinant_df[r1_recombinant_df.loc[:, 'big_fail'] == False]
r2_true_recombinants_df = r2_recombinant_df[r2_recombinant_df.loc[:, 'big_fail'] == False]
# r1_recombinant_df[r1_recombinant_df['in_guide_library'] == True]

# print(len(r1_true_recombinants_df))
# print(len(r2_true_recombinants_df))


big_fail_pcnt = round((N_big_fails/recombinant)*100, 2)

num_r1_failed_recombinants = r1_failed_recombinants['read_group']
num_r2_failed_recombinants = r2_failed_recombinants['read_group']

print("#"*46)
print(f"Of the {recombinant} recombinant read groups, {N_big_fails} read groups had a sequence not in the guide list, so {big_fail_pcnt} % of recombinants can be considered failures.")
print("#"*46)

r1_failed_recombinants['read_group']

# len(r1_failed_recombinants['read_group'])

# Start stacking the recombinants!
r1_failed_recombinant_stacked_df = r1_failed_recombinants_df[['title', 'seq', 'plus', 'qual']].stack()
r2_failed_recombinant_stacked_df = r2_failed_recombinants_df[['title', 'seq', 'plus', 'qual']].stack()

r1_failed_recombinant_fastq_df = r1_failed_recombinant_stacked_df
r2_failed_recombinant_fastq_df = r2_failed_recombinant_stacked_df

r1_true_recombinant_stacked_df = r1_true_recombinants_df[['title', 'seq', 'plus', 'qual']].stack()
r2_true_recombinant_stacked_df = r2_true_recombinants_df[['title', 'seq', 'plus', 'qual']].stack()

r1_true_recombinant_fastq_df = r1_true_recombinant_stacked_df
r2_true_recombinant_fastq_df = r2_true_recombinant_stacked_df

# Export the new fastq.
# Define the child directory relative to the current script's location. Use lines for shell script
# script_dir = os.path.dirname(os.path.abspath(__file__))  # Absolute path to the script
# child_directory = os.path.join(script_dir, 'output')

# Define the child directory relative to the current working directory
child_directory = os.path.join(os.getcwd(), 'output')

# Ensure the child directory exists
if not os.path.exists(child_directory):
    os.makedirs(child_directory)

# Extract just the filename from r1_file and r2_file
r1_filename = os.path.basename(r1_file)
r2_filename = os.path.basename(r2_file)

print(f"Saving files to directory: {child_directory}")

# Now set the output file paths in the subdirectory 'hits'
r1_hits_out_file = os.path.join(child_directory, "hits." + r1_filename)
r2_hits_out_file = os.path.join(child_directory, "hits." + r2_filename)
r1_true_recombinant_out_file = os.path.join(child_directory, "recombinants." + r1_filename)
r2_true_recombinant_out_file = os.path.join(child_directory, "recombinants." + r2_filename)
r1_failed_recombinant_out_file = os.path.join(child_directory, "fails." + r1_filename)
r2_failed_recombinant_out_file = os.path.join(child_directory, "fails." + r2_filename)

# Save the DataFrames to their respective CSV files
r1_hits_fastq_df.to_csv(r1_hits_out_file, sep='\t', index=False, header=False, compression='gzip')
r2_hits_fastq_df.to_csv(r2_hits_out_file, sep='\t', index=False, header=False, compression='gzip')
r1_true_recombinant_fastq_df.to_csv(r1_true_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
r2_true_recombinant_fastq_df.to_csv(r2_true_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
r1_failed_recombinant_fastq_df.to_csv(r1_failed_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
r2_failed_recombinant_fastq_df.to_csv(r2_failed_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')



# # Define the child directory relative to the current script's location
# child_directory = os.path.join(os.path.dirname(__name__), 'output')
#
# # Ensure the child directory exists
# if not os.path.exists(child_directory):
#     os.makedirs(child_directory)
#
# # Set the output file paths in the child directory
# r1_hits_out_file = os.path.join(child_directory, "hits." + r1_file)
# r2_hits_out_file = os.path.join(child_directory, "hits." + r2_file)
# r1_true_recombinant_out_file = os.path.join(child_directory, "recombinants." + r1_file)
# r2_true_recombinant_out_file = os.path.join(child_directory, "recombinants." + r2_file)
# r1_failed_recombinant_out_file = os.path.join(child_directory, "fails." + r1_file)
# r2_failed_recombinant_out_file = os.path.join(child_directory, "fails." + r2_file)
#
# # Save the DataFrames to their respective CSV files
# r1_hits_fastq_df.to_csv(r1_hits_out_file, sep='\t', index=False, header=False, compression='gzip')
# r2_hits_fastq_df.to_csv(r2_hits_out_file, sep='\t', index=False, header=False, compression='gzip')
# r1_true_recombinant_fastq_df.to_csv(r1_true_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
# r2_true_recombinant_fastq_df.to_csv(r2_true_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
# r1_failed_recombinant_fastq_df.to_csv(r1_failed_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')
# r2_failed_recombinant_fastq_df.to_csv(r2_failed_recombinant_out_file, sep='\t', index=False, header=False, compression='gzip')

"""# Add some narrative.

"""

print("#"*46)
print("Your analysis has just finished.")
print("Reads from matched read groups on whose guides were on target for both R1 and R2 are found in the files prefixed hits.*.")
print("Reads from matched read groups on whose guides were on were off target and considered recombinant for either R1 and R2 are found in the files prefixed recombinants.*.")
print("Reads from recombinant read groups on whose guides were not mathced to a known guide sequence for either R1 and R2 are found in the files prefixed fails.*.")
print("Good luck and feel free to generally ignore any outputs below here!")
print("#"*46)

"""# Just for executable testing below."""

# Commented out IPython magic to ensure Python compatibility.
# %%bash
# # python 19bp_dual_guide_parser_tool.py --help
# python 19_20_var_bp_dual_guide_parser.py --guides_file 20200513_library_1_2_unbalanced_dJR051.csv --r1_file UDP0007_S1_R1_001.fastq.gz --r2_file UDP0007_S1_R2_001.fastq.gz --N_reads 100000
