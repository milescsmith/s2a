from collections import OrderedDict
from typing import List, Optional

import numpy as np
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse.csc import csc_matrix


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
    adata: AnnData, var_features: List[str], meta_features: DataFrame = None
) -> AnnData:

    if isinstance(adata.X, csc_matrix):
        exprs = adata.X.toarray()
    elif isinstance(adata.X, np.ndarray):
        exprs = adata.X

    if meta_features is not None and len(meta_features.columns) > 0:
        if "sct.gmean" in meta_features.columns:
            var_dict = {
                "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, exprs),
                "highly_variable": [
                    True if _ in var_features else False for _ in adata.var.index
                ],
                "means": meta_features["sct.gmean"],
                "dispersions": meta_features["sct.variance"]
                if "sct.variance" in meta_features.columns
                else np.zeros(len(meta_features)),
                "residuals_dispersions": meta_features["sct.residual_variance"]
                if "sct.residual_variance" in meta_features.columns
                else np.zeros(len(meta_features)),
                "residual_mean": meta_features["sct.residual_mean"]
                if "sct.residual_mean" in meta_features.columns
                else np.zeros(len(meta_features)),
                "detection_rate": meta_features["sct.detection_rate"]
                if "sct.detection_rate" in meta_features.columns
                else np.zeros(len(meta_features)),
            }
        elif "vst.mean" in meta_features.columns:
            var_dict = {
                "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, exprs),
                "highly_variable": [
                    True if _ in var_features else False for _ in adata.var.index
                ],
                "means": meta_features["vst.mean"],
                "dispersions": meta_features["vst.variance"]
                if "vst.variance" in meta_features.columns
                else np.zeros(len(meta_features)),
                "dispersions_norm": meta_features["vst.variance.standardized"]
                if "vst.variance.standardized" in meta_features.columns
                else np.zeros(len(meta_features)),
            }
        else:
            var_dict = {
                "n_cells": np.apply_along_axis(lambda x: (x > 0).sum(), 0, exprs),
                "highly_variable": [
                    True if _ in var_features else False for _ in adata.var.index
                ],
                "means": meta_features["mean"]
                if "mean" in meta_features.columns
                else np.zeros(len(meta_features)),
                "dispersions": meta_features["variance"]
                if "variance" in meta_features.columns
                else np.zeros(len(meta_features)),
                "dispersions_norm": meta_features["variance.standardized"]
                if "variance.standardized" in meta_features.columns
                else np.zeros(len(meta_features)),
            }
        adata.var = DataFrame.from_dict(data=var_dict)

    return adata


def add_reduction(
    adata: AnnData,
    reduction_key: str,
    cell_embeddings: np.ndarray,
    feature_loadings: Optional[np.ndarray] = None,
    reduction_sd: Optional[np.ndarray] = None,
) -> AnnData:

    adata.obsm[f"X_{reduction_key}"] = cell_embeddings
    if feature_loadings is not None and feature_loadings.size > 0:
        adata.varm["PCs"] = feature_loadings

    if reduction_sd is not None and len(reduction_sd) > 0:
        reduc = OrderedDict()
        reduc["variance"] = reduction_sd
        adata.uns[reduction_key] = reduc

    return adata
