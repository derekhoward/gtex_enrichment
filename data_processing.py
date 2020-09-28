import numpy as np
import pandas as pd
from scipy.stats import zscore

print('Loading expression data...')
metadf = pd.read_table(
    './data/raw/GTEx_v1.9/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')
exp = pd.read_csv('./data/processed/gtex_sample_expression.csv',
                  index_col='gene_symbol')

exp = exp.applymap(lambda x: np.log(x+1))


# extract sample metadata
sample_metadata = metadf[metadf.SAMPID.isin(exp.columns)].loc[:, [
    'SAMPID', 'SMTS', 'SMTSD', 'SMUBRID']]
sample_metadata['donor_id'] = sample_metadata.SAMPID.str.split(
    '-').str[:2].str.join('-')
print('Merge expression with metadata')
merged_data = exp.T.merge(sample_metadata.loc[:, [
    'SAMPID', 'donor_id', 'SMTSD']], left_index=True, right_on='SAMPID')
merged_data = merged_data.set_index('SAMPID')

donors = merged_data.donor_id.unique()
print(f'There are {len(donors)} donors in the dataset')


# SMTSD is tissue sampled
def process_donor_subset(subset_df):
    meta = subset_df.loc[:, ['donor_id', 'SMTSD']]
    ranked_genes_in_sample = subset_df.drop(
        ['donor_id', 'SMTSD'], axis=1).rank(axis=1)
    zscored_across_samples = ranked_genes_in_sample.apply(zscore)
    processed = pd.concat([zscored_across_samples, meta], axis=1)
    return processed


print('Processing donors individually')
transformed_donor_dfs = []
for donor in donors:
    subset = merged_data[merged_data.donor_id == donor]
    transformed_subset = process_donor_subset(subset)
    transformed_donor_dfs.append(transformed_subset)

all_donors = pd.concat(transformed_donor_dfs)
tissue_level_mean = all_donors.drop('donor_id', axis=1).groupby('SMTSD').mean()
tissue_level_mean.to_csv('./data/processed/gtex_processed.csv')
tissue_level_ranks = tissue_level_mean.T.rank()
tissue_level_ranks.index.name = 'gene_symbol'
tissue_level_ranks.to_csv('./data/processed/gtex_processed_ranks.csv')
