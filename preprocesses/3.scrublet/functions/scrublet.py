import sys
from pandas import Index
import scanpy as sc
from scrublet import Scrublet
import numpy as np

sample_name = sys.argv[1]
adata = sc.read_10x_h5(
    sys.argv[2]
)
adata.var_names_make_unique()
adata.obs_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=5)
df = adata.to_df().astype(int)
print(df.index)
print(df.columns)
scrub = Scrublet(df)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
print(doublet_scores)
len_all = 0
if predicted_doublets is not None:
    len_all = len(df[predicted_doublets].index)
    df = df[predicted_doublets == False]
print(f"{len_all} cells were removed.")
df.index = Index(
    np.array(
        [f"{sample_name}_{str(x)}" for x in range(
            len(df.index))]).tolist(),
    dtype=str
)
print(df.index)
print(df.columns)
df.T.to_csv(
    sys.argv[3],
    sep="\t",
    header=True
)
