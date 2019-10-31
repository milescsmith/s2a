from anndata import AnnData
from pandas import DataFrame
import numpy as np
from typing import List, Optional
from collections import OrderedDict


def add_meta_data(
    adata: AnnData, md: DataFrame, column_for_rownames: Optional[str] = None
) -> AnnData:
    if "Unnamed: 0" in md.columns:
        md.index = md["Unnamed: 0"]
        md.drop(["Unnamed: 0"], axis=1)
    if column_for_rownames:
        md.index = md[column_for_rownames]
        md.drop([column_for_rownames], axis=1)
    adata.obs = md
    return adata


def add_feature_data(
    adata: AnnData, var_features: List[str], meta_features: DataFrame
) -> AnnData:

    if "sct.gmean" in meta_features.columns:
        var_dict = {
            "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, adata.X),
            "highly_variable": [
                True if _ in var_features else False for _ in adata.var.index
            ],
            "means": meta_features["sct.gmean"],
            "dispersions": meta_features["sct.variance"],
            "residuals_dispersions": meta_features["sct.residual_variance"],
            "residual_mean": meta_features["sct.residual_mean"],
            "detection_rate": meta_features["sct.detection_rate"],
        }
    elif "vst.mean" in meta_features.columns:
        var_dict = {
            "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, adata.X),
            "highly_variable": [
                True if _ in var_features else False for _ in adata.var.index
            ],
            "means": meta_features["vst.mean"],
            "dispersions": meta_features["vst.variance"],
            "dispersions_norm": meta_features["vst.variance.standardized"],
        }
    else:
        var_dict = {
            "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, adata.X),
            "highly_variable": [
                True if _ in var_features else False for _ in adata.var.index
            ],
            "means": meta_features["mean"],
            "dispersions": meta_features["variance"],
            "dispersions_norm": meta_features["variance.standardized"],
        }
    adata.var = DataFrame.from_dict(data=var_dict)
    return adata


def add_reduction(
    adata: AnnData,
    reduction_key: str,
    cell_embeddings: np.ndarray,
    feature_loadings: np.ndarray,
    reduction_sd: np.ndarray,
) -> AnnData:

    adata.obsm[f"X_{reduction_key}"] = cell_embeddings
    if feature_loadings.size > 0:
        adata.varm["PCs"] = feature_loadings
    if len(reduction_sd) > 0:
        reduc = OrderedDict()
        reduc["variance"] = reduction_sd
        adata.uns[reduction_key] = reduc
    return adata
