# CSCC-SingleCell
  
### This directory contains the final files used for integration and analysis 

#### File Contents::

AUCell_score_MDK_up_down_genes_sig.Rmd - R notebook used to calculate the probability of cells being TSK using AUCell and genes from Ji et al. This was unsucessful.

Cell Type annotation Full Resolution mt-check.ipynb - Jupyter notebook Containing downstream analysis, annotations, and UMAP plots of the final anndata object

Celltypist annotations.ipynb - Testing Celltypist package on each of the Zou, Ji, Yost, and Lyko Datasets individually

DE_prep_for_GSEA.ipynb - Running sc.tl.rank_gene_groups on different keratinocyte subsets to produce ranked gene files for use in GSEA  analysis

Full_object_final.ipynb - Running my full pre-processing and integration pipeline on the combined object, using the script SCProcess.py

Individual-Tumor.ipynb - Running SCProcess.py on all 4 datasets individually to visualize clustering before combining them

Ji, Lyko, Yost, Zou files -- Combining the diles from GEO with the metadata and transposing to .h5ad files. Standardizing data to include batch, patient and Condition columns

SCProcess.py - Pre-processing and integration pipeline with parameters used in the final object

