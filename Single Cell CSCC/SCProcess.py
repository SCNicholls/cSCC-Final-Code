import scanpy as sc
import harmonypy as hm
import pandas as pd
import anndata as ad
import numpy as np
import scrublet as scr
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

random_State=7

# %%
def QC(adata):
    #Doing quality control on the data to identify filtering thresholds
    # find genes that are mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # find genes that are ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=None
    )
    print(adata.obs[
        ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", "pct_counts_ribo"]
    ].describe())
    print(adata.var.describe())
    
    # violin plot of some distributions
    sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.2,
        multi_panel=True,
        rotation=30,
        save=False,
        show=True,
    )
    #### Filter cells according to identified QC thresholds:
    print("Total number of cells: {:d}".format(adata.n_obs))

    # remove cells with more than 20% MT genes
    adata = adata[adata.obs.pct_counts_mt < 30, :].copy()
    print("Number of cells after mt filter: {:d}".format(adata.n_obs))
    sc.pp.filter_cells(adata, min_genes=300)
    print("Number of cells after gene filter: {:d}".format(adata.n_obs))

    #Removing mitochondrial and ribosomal genes
    print("Total number of genes: {:d}".format(adata.shape[1]))
    adata = adata[:, ~adata.var["mt"].values]
    print("After mt genes removal: ", adata.shape[1])
    adata = adata[:, ~adata.var["ribo"].values]
    print("After ribo genes removal: ", adata.shape[1])

    adata.write_h5ad("/data/scratch/ha20577/adata_after_QC.h5ad")
    return adata

# %%
def DoubletRM(adata):
    #Removing doublets
    #Using Scrublet to identify doublets
    #Do I need to use batch aware processing?
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["doublet_scores"] = doublet_scores
    adata.obs["predicted_doublets"] = predicted_doublets

    #Filtering out the doublets
    adata = adata[~adata.obs["predicted_doublets"], :].copy()
    print("After removing doublets: ", adata.shape[0], "cells")

    adata.write_h5ad("/data/scratch/ha20577/adata_after_doubletRM.h5ad")

    return adata

 

# %%
def Normalize(adata):
    # Normalizing the data and keeping it raw
    adata.layers["raw"] = adata.X.copy()  # preserve counts

    # normalize + log1p
    # Total-count normalize (library-size correct) the data matrix
    # to 10,000 reads per cell, so that counts become comparable among cells.
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    adata.layers["normalised"] = adata.X.copy()
    # Logarithmize the data
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    adata.raw = adata.copy()  # keep normalised log1p
    return adata
    
# %%
def Reduction(adata):

    sc.pp.highly_variable_genes(
        adata,
        subset=True,  # subset for integration (but full lognorm data in .raw)
        layer="log1p",
        flavor="seurat_v3",
        n_top_genes=2000,
        span=0.3,
        min_disp=0.5,
        min_mean=0.0125,
        max_mean=3
        # no batch correction as we are integrating
    )


#sc.pp.highly_variable_genes
    
    print("AnnData HVG dimensions: ", adata.X.shape)
    print("AnnData Raw dimensions: ", adata.raw.X.shape)

    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

    # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)

    #PCA
    sc.tl.pca(adata, svd_solver='arpack', random_state=random_State)
    ### Scatter plot for PCA, but we will not use later on
    sc.pl.pca(adata, color=['n_genes_by_counts','pct_counts_mt'], color_map='viridis')
    ### Estimate number of PCs to use: (rough estimate is often fine)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

    adata.write_h5ad("/data/scratch/ha20577/adata_after_reduction.h5ad") 

# %%
def Integration(adata, cat_var_to_regress='batch'):
    #Integrating the data using Harmony
    #Run Harmony
    #cat_var_to_regress must be string
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, cat_var_to_regress, max_iter_harmony=20)
    adata.obsm['X_pca_harmony'] = np.transpose(ho.Z_corr) # add corrected Harmony principal components to 
    adata.write_h5ad("/data/scratch/ha20577/adata_after_integration.h5ad") 


def Neighborhood(adata):
    random_State = 7
    ### Construct neighbourhood graph using corrected principal components generated with Harmony
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep='X_pca_harmony', random_state=random_State)
    ### embed neighbourhood graph
    sc.tl.umap(adata, random_state=random_State)

    ### Plot umaps of quality control metrics
    sc.pl.umap(adata, color=['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo'], color_map='viridis',
            vmax = [15000,6000,20,50,0.5], ncols=5)
    ### 
    sc.pl.umap(adata, color=["Condition", 'Patient_ID'], color_map="viridis")

    #Clustering to determine resolution
    sc.tl.leiden(adata, resolution=0.6, key_added = 'clusters_r06', random_state=random_State)
    sc.tl.leiden(adata, resolution=0.8, key_added = 'clusters_r08', random_state=random_State)
    sc.tl.leiden(adata, resolution=1, key_added = 'clusters_r1', random_state=random_State)

    sc.pl.umap(adata, color = ['clusters_r06','clusters_r08','clusters_r1'])

    #Saving the file thus far in case it crashes during DE analysis
    #adata.write_h5ad("/data/scratch/ha20577/Combined_adata.h5ad") 

    adata.write_h5ad("/data/scratch/ha20577/adata_after_neighborhood.h5ad") 

# %%
def DE(adata, cluster_key='clusters_r08'):
    #Differential expression analysis
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon', pts=True, use_raw=True)
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
    
    pval_thresh = 0.05
    log2fc_thresh = 0.25
    pct_cutoff = 0.1
    cluster_de_genes = dict()
    for cluster in sorted(set(adata.obs[cluster_key])):
        #Get DE genes for each cluster
        cluster_de_genes[cluster] = sc.get.rank_genes_groups_df(adata,
                                                                group=cluster, 
                                                                key='rank_genes_groups', 
                                                                pval_cutoff=pval_thresh, 
                                                                log2fc_min=log2fc_thresh, 
                                                                log2fc_max=None).sort_values('logfoldchanges',ascending=False)
        cluster_de_genes[cluster] = cluster_de_genes[cluster][cluster_de_genes[cluster]['pct_nz_group'] > pct_cutoff]


# %%
def Runall(adata):
    adata = QC(adata)
    adata = DoubletRM(adata)
    adata = Normalize(adata)
    Reduction(adata)
    Integration(adata)
    Neighborhood(adata)
    DE(adata)
    adata.write_h5ad("/data/BCI-SingleCell/SCC_Atlas/Sam_Nicholls/Combined_adata_object_HR.h5ad") 
    return adata


