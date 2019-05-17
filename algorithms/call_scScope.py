'''The goal of this demo is to show how to identify cell subpopulations based on latent
representations of gene expression learned by scScope.'''
import scscope_small as DeepImpute
import pandas as pd
import phenograph
import pickle
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

# For this demo we normalize data using scanpy which is not a required package for scScope.
# To install, use: pip install scanpy
import scanpy.api as sc


if __name__ == '__main__':
    import time
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts_file', metavar='counts_file', default=None, required=True, help='counts_file file')
    parser.add_argument('--num_pc', metavar='N', type=int, default=None, help='num_pc file')
    parser.add_argument('--out', metavar='out', required=True, help='output file')
    args = parser.parse_args()

    # 1. Load gene expression matrix of simulated data
    counts_drop = pd.read_csv(args.counts_file, header=0, index_col=0)
    # 2. Normalize gene expression based on scanpy (normalize each cell to have same library size)
    # matrix of cells x genes
    start_time = time.time()
    gene_expression = sc.AnnData(counts_drop.values)
    # normalize each cell to have same count number
    sc.pp.normalize_per_cell(gene_expression)
    # update datastructure to use normalized data
    gene_expression = gene_expression.X

	# 3. scScope learning
    if gene_expression.shape[0] >= 100000:
        DI_model = DeepImpute.train(gene_expression, args.num_pc, T=2, batch_size=512, max_epoch=10, num_gpus=4)
    else:
        DI_model = DeepImpute.train(gene_expression, args.num_pc, T=2, batch_size=64, max_epoch=300, num_gpus=4)

    # 4. latent representations and imputed expressions
    latent_code, imputed_val, _ = DeepImpute.predict(gene_expression, DI_model)
    elapsed_time = time.time() - start_time

    np.savetxt(args.out, latent_code, fmt='%.4f')
