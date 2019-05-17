s2a <- NULL
anndata <- NULL

.onLoad <- function(libname, pkgname) {
  s2a <<- reticulate::import(module = "s2a", delay_load = TRUE)
  anndata <<- reticulate::import(module = "anndata", delay_load = TRUE)
  `%nin%` <- purrr::compose(`!`, `%in%`)
}
