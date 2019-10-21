library(reticulate)
use_condaenv("reticulate")
load("data/pbmc3k.RData")

test_that("conversion works", {
  pbmc3k_adata = convert_to_anndata(object = pbmc3k, assay = "RNA")
  expect_equal(dim(pbmc3k_adata$X), c(2700, 32738))
  expect_equal(dim(pbmc3k_adata$obsm["X_pca"]), c(2700, 50))
})
