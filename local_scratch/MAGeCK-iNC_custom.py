import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import random
from random import shuffle
import subprocess
import shlex
import sys
import pandas as pd
from debugpy._vendored.pydevd.pydev_sitecustomize.sitecustomize import input

pd.set_option('display.max_rows', 150)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1000)

def product_threshold_fdr(df, fdr = 0.05):
    maxi = abs(df['product']).max()
    for pro in np.arange(0, maxi, 0.1):
        df_thres = df[abs(df['product']) > pro]
        if (1.0 * len(df_thres[df_thres['index'].str.contains('negative_control')]) / len(df_thres)) < fdr:
            break
    return pro, df_thres


def rank_test(df):
    df = df[df['treat_mean'] > 20]
    df_ntc = df[df['Gene'].str.contains('negative_control')]
    df_targeting = df[~df['Gene'].str.contains('negative_control')]
    df_targeting = df_targeting[df_targeting.duplicated('Gene', keep=False)]
    ntc_sgRNA_p = list(df_ntc['p.twosided'])
    ntc_sgRNA_p_lfc = zip(list(df_ntc['p.twosided']), list(df_ntc['LFC']))
    genes = df_targeting['Gene'].unique()
    num_of_genes = len(genes)
    gene_lfc_p = {}

    for gene in genes:
        df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
        lfc = df_gene.iloc[:3]['LFC'].mean()
        x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3], ntc_sgRNA_p, alternative = 'two-sided')
        gene_lfc_p[gene] = [lfc, pvalue]

    random.seed(10)
    for j in range(num_of_genes):
        ntc_sgRNA_p_lfc = list(ntc_sgRNA_p_lfc)
        shuffle(ntc_sgRNA_p_lfc)
        ntc_selected = ntc_sgRNA_p_lfc[:5]
        ntc_selected_p = [i[0] for i in ntc_selected]
        ntc_lfc = np.mean([i[1] for i in sorted(ntc_selected, key = lambda x: x[0])][:3])
        x, ntc_pvalue = mannwhitneyu(ntc_selected_p, ntc_sgRNA_p, alternative = 'two-sided')
        gene_lfc_p['negative_control_' + str(j)] = [ntc_lfc, ntc_pvalue]
    return gene_lfc_p

# u test
# df_mageck = pd.read_table(output_folder + "/" + output_name + '.sgrna_summary.txt')
# df_mageck = pd.read_table('/Users/Claire/Downloads/thresh_300_.sgrna_summary.txt')
df_mageck = pd.read_table('/Users/Claire/Library/CloudStorage/Box-Box/Graduate_School/PhD/00_projects/CRISPRi_alg_optimization/00_references_and_code/files_from_others/01_NINDS_informatics/3_mageck_vispr_analysis/test/High_vs_Low.sgrna_summary.txt')
df = pd.DataFrame(rank_test(df_mageck)).T
df_copy = df.copy()
df.columns = ['epsilon', 'pvalue']
df.reset_index(inplace = True)
df['gene'] = df['index'].apply(lambda x: x.split('_')[0])
df['epsilon'] = -df['epsilon']
df_ntc = df[df['gene'] == 'negative']
df['sub_ntc_med_epsilon'] = df['epsilon'] - df_ntc['epsilon'].median()
# df['epsilon'] = df['epsilon'] - df_ntc['epsilon'].median()
df['product'] = df['sub_ntc_med_epsilon'] * (-np.log10(df['pvalue']))
df.sort_values('product', ascending = False)

# FDR
thres, df_hits = product_threshold_fdr(df, fdr=0.1)
# df.sort_values('product', ascending = False).to_csv(output_folder + '/' + output_name + '_all_genes.csv', index = False)
# df_hits.sort_values('product', ascending = False).to_csv(output_folder + '/' + output_name + '_fdr%s_product%s_hits.csv' % (fdr, thres),
#                                                          index = False)
df_ntc = df[df['index'].str.contains('negative_control')]


npg = ["#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"]
plt.figure(figsize = [10, 8])
df_pos = df_hits[df_hits['epsilon'] > 0]
df_neg = df_hits[df_hits['epsilon'] < 0]
plt.scatter(df_pos['epsilon'], -np.log10(df_pos['pvalue']), c = "#DC0000FF", s = 5, label = 'Positive hits')
plt.scatter(df_neg['epsilon'], -np.log10(df_neg['pvalue']), c = "#3C5488FF", s = 5, label = 'Negative hits')
plt.scatter(df['epsilon'], -np.log10(df['pvalue']), c = '#F39B7FFF', s = 5, label = 'Other genes')
plt.scatter(df_pos['epsilon'], -np.log10(df_pos['pvalue']), c = "#DC0000FF", s = 5, label = None)
plt.scatter(df_neg['epsilon'], -np.log10(df_neg['pvalue']), c = "#3C5488FF", s = 5, label = None)
plt.scatter(df_ntc['epsilon'], -np.log10(df_ntc['pvalue']), c = 'grey', s = 5, label = "Negative control")
genes = list(df['gene'])
phenotype = list(df['epsilon'])
p_value = list(-np.log10(df['pvalue']))
i = 0
# for x, y, s in zip(phenotype, p_value, genes):
#     if s in genes_to_label:
#         plt.annotate(s, (x, y), fontsize = 16)
#         if i == 0:
#             plt.scatter(x, y, c = 'darkgreen', s = 20, label = 'Genes of interest')
#         else:
#             plt.scatter(x, y, c = 'darkgreen', s = 20)
#         i = 1

    # genes = list(df_hits['gene'])
    # phenotype = list(df_hits['epsilon'])
    # p_value = list(-np.log10(df_hits['pvalue']))
    # for x, y, s in zip(phenotype, p_value, genes):
    #     plt.annotate(s, (x, y), fontsize = 10)
x = np.arange(0.01, 10, 0.01)
y = [thres / i for i in x]
plt.plot(x, y, '--', c = 'k', )
plt.plot(-x, y, '--', c = 'k', label = 'FDR = %s' % 0.1)
lim = max(abs(df['epsilon'].min() - 1), (df['epsilon'].max() + 2))
plt.xlim(-lim, lim)
plt.ylim(0, -np.log10(df['pvalue'].min()) + 0.5)
plt.legend(loc = 1, fontsize = 'large', fancybox = True)
plt.xlabel('Phenotype', fontsize = 14)
plt.ylabel('-log10 P', fontsize = 14)
# plt.title(output_name, fontsize = 18)
# plt.savefig(output_folder + "/" + output_name + '_volcano_plot.pdf')
# Label the points with gene names
for i, row in df_pos.iterrows():
    plt.annotate(row['gene'], (row['epsilon'], -np.log10(row['pvalue'])), fontsize=8, color='black', alpha=0.8)

for i, row in df_neg.iterrows():
    plt.annotate(row['gene'], (row['epsilon'], -np.log10(row['pvalue'])), fontsize=8, color='black', alpha=0.8)
plt.show()






import matplotlib.pyplot as plt
import numpy as np
import mplcursors

# Assuming df, df_hits, and df_ntc are already defined dataframes
npg = ["#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"]

plt.figure(figsize=[10, 8])

df_pos = df_hits[df_hits['epsilon'] > 0]
df_neg = df_hits[df_hits['epsilon'] < 0]

# Scatter plots for different categories
scatter_pos = plt.scatter(df_pos['epsilon'], -np.log10(df_pos['pvalue']), c="#DC0000FF", s=5, label='Positive hits')
scatter_neg = plt.scatter(df_neg['epsilon'], -np.log10(df_neg['pvalue']), c="#3C5488FF", s=5, label='Negative hits')
scatter_other = plt.scatter(df['epsilon'], -np.log10(df['pvalue']), c='#F39B7FFF', s=5, label='Other genes')
scatter_ntc = plt.scatter(df_ntc['epsilon'], -np.log10(df_ntc['pvalue']), c='grey', s=5, label='Negative control')

# FDR threshold lines
x = np.arange(0.01, 10, 0.01)
y = [thres / i for i in x]
plt.plot(x, y, '--', c='k')
plt.plot(-x, y, '--', c='k', label='FDR = 0.1')

# Set plot limits
lim = max(abs(df['epsilon'].min() - 1), (df['epsilon'].max() + 2))
plt.xlim(-lim, lim)
plt.ylim(0, -np.log10(df['pvalue'].min()) + 0.5)

# Add labels and legend
plt.xlabel('Phenotype', fontsize=14)
plt.ylabel('-log10 P', fontsize=14)
plt.legend(loc=1, fontsize='large', fancybox=True)

# Use mplcursors to add interactive annotations
mplcursors.cursor([scatter_pos, scatter_neg, scatter_other, scatter_ntc], hover=True).connect(
    "add", lambda sel: sel.annotation.set_text(f'{df["gene"].iloc[sel.target.index]}'))

# Show plot
plt.show()


