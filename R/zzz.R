s2a <- NULL
anndata <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  py_config <- reticulate::py_discover_config(required_module = "s2a")
  s2a <- reticulate::import(module = "s2a", delay_load = TRUE)
  anndata <- reticulate::import(module = "anndata", delay_load = TRUE)
  `%nin%` <- purrr::compose(`!`, `%in%`)
}