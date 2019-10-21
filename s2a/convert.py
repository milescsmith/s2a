from anndata import AnnData
from pandas import DataFrame
import numpy as np
from typing import List


def add_meta_data(adata: AnnData, 
                  md: DataFrame) -> AnnData:
    md.index = md["Unnamed: 0"]
    md.drop(["Unnamed: 0"], axis=1)
    adata.obs = md
    return adata


def add_feature_data(adata: AnnData,
                     var_features: List[str],
                     meta_features: DataFrame) -> AnnData:
    
    var_dict = {"n_cells": np.apply_along_axis(lambda x: (x>0).sum(), 1, adata.X),
                "highly_variable": [True if _ in var_features else False for _ in adata.var.index],
                "means": meta_features["means"],
                "dispersions": meta_features['variance'],
                "dispersions_norm": meta_features["variance.standardized"]}
    adata.var = DataFrame.from_dict(data = var_dict)
    return adata


def add_reduction(adata: AnnData,
                   reduction_key: str,
                   cell_embeddings: np.ndarray,
                   feature_loadings: np.ndarray,
                   reduction_sd: np.ndarray) -> AnnData:
    
    adata.obsm[f"X_{reduction_key}"] = cell_embeddings
    if feature_loadings.size > 0:
        adata.varm["PCs"] = feature_loadings
    if reduction_sd.size > 0:
        adata.uns[reduction_key]["variance"] = reduction_sd
    return adata