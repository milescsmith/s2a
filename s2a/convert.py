from anndata import AnnData
from pandas import DataFrame
from scanpy import tl
import numpy as np
from typing import List


def add_meta_data(adata: AnnData, 
                  md: DataFrame) -> AnnData:
    adata.obs = md
    if "seurat_clusters" in md.keys():
        adata.obs["louvain"] = md["seurat_clusters"]
    else:
        tl.louvain(adata)
    return adata


def add_feature_data(adata: AnnData,
                     var_features: List[str],
                     meta_features: DataFrame) -> AnnData:
    if "gmean" in meta_features.keys():
        var_dict = {"n_cells": np.apply_along_axis(lambda x: (x>0).sum(), 0, adata.X),
                    "highly_variable": [True if _ in var_features else False for _ in adata.var.index],
                    "means": meta_features["gmean"],
                    "dispersions": meta_features['variance'],
                    "dispersions_norm": meta_features["residual_variance"],
                    "detection_rate": meta_features["detection_rate"]}
    else:
        var_dict = {"n_cells": np.apply_along_axis(lambda x: (x>0).sum(), 0, adata.X),
                    "highly_variable": [True if _ in var_features else False for _ in adata.var.index],
                    "means": meta_features["mean"],
                    "dispersions": meta_features['variance'],
                    "dispersions_norm": meta_features["variance.standardized"]}
    adata.var = DataFrame.from_dict(data = var_dict)
    return adata


def add_reduction(adata: AnnData,
                   reduction_key: str,
                   cell_embeddings: np.ndarray,
                   feature_loadings: np.ndarray,
                   reduction_sd: np.array) -> AnnData:
    
    adata.obsm[f"X_{reduction_key}"] = cell_embeddings
    if np.size(feature_loadings) > 0:
        adata.varm["PCs"] = feature_loadings
    if len(reduction_sd) > 0:
        variance_dict = {"variance" : np.array(reduction_sd),
                         "variance_ratio" : np.array(reduction_sd/reduction_sd.sum())}
        adata.uns[reduction_key] = variance_dict
    return adata

