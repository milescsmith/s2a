#' @title s2a
#' @description Convert Seurat object to AnnData object
#'
#' @param object Seurat object to convert
#' @param assay Seurat object assay to add to new object. Default: "RNA"
#' @param slot Seurat assay slot to convert. If `data` or `scale.data` is used, 
#' the `counts` slot (if available) is added to `adata.raw`
#' 
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom Matrix t
#' @importFrom purrr walk
#' @importFrom reticulate import
#' @importFrom rlang %||%
#' @importFrom Seurat GetAssayData
#' @importFrom tibble tibble as_tibble column_to_rownames
#'
#' @return anndata$AnnData object
#' @export
#' 
convert_to_anndata <-
  function(
    object,
    assay = NULL,
    slot  = "data"
    ){
    
  s2a <-
    import(
      module     = "s2a",
      delay_load = TRUE
      )
  
  anndata <-
    import(
      module     = "anndata",
      delay_load = TRUE
      )
  
  adata <- anndata$AnnData()
  
  assay <- assay %||% DefaultAssay(object)
  
  exprDat <-
    GetAssayData(
      object = object, 
      assay  = assay,
      slot   = "data"
      )
  
  texprDat <- Matrix::t(exprDat)
  
  adata <- anndata$AnnData(texprDat)

  adata$obs_names <- rownames(texprDat)
  adata$var_names <- colnames(texprDat)
  
  if (slot != "counts" & "counts" %in% names(object[[assay]])) {
    raw <-
      GetAssayData(
        object = object,
        assay  = assay,
        slot   = "counts"
        )
    traw      <- t(raw)
    adata$raw <- traw
  }

  adata <-
    s2a$add_meta_data(
      adata = adata, 
      md    = object@meta.data
      )
  
  if (ncol(object[[assay]]@meta.features) == 0){
    adata <-
      s2a$add_feature_data(
        adata         = adata,
        var_features  = object[[assay]]@var.features,
        meta_features = NULL
      )
  } else {
    adata <-
      s2a$add_feature_data(
        adata         = adata,
        var_features  = object[[assay]]@var.features,
        meta_features = object[[assay]]@meta.features
      )
  }
  
  
  for (i in names(object@reductions)){
    print(i)
    if ((nrow(object[[i]]@feature.loadings) > 0           ) & 
        (nrow(object[[i]]@feature.loadings) < nrow(object))){
      
      
      feature_loadings = object[[i]]@feature.loadings
      new_loadings <-
        tibble::tibble(
          feature = rownames(object)[!rownames(object) %in% rownames(feature_loadings)]
          )
      
      new_loadings[,colnames(feature_loadings)] <- 0
      
      feature_loadings <-
        feature_loadings |>
        tibble::as_tibble(rownames = "feature") |>
        rbind(new_loadings) |>
        dplyr::arrange(feature) |>
        tibble::column_to_rownames('feature') |>
        as.matrix()
    } else {
      feature_loadings <- object[[i]]@feature.loadings
    }
    try(
      adata <-
        s2a$add_reduction(
          adata            = adata,
          reduction_key    = i,
          cell_embeddings  = object[[i]]@cell.embeddings,
          feature_loadings = feature_loadings,
          reduction_sd     = object[[i]]@stdev
          )
      )
  }
  
  return(adata)
}
