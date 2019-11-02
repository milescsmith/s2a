.onLoad <- function(libname, pkgname) {
  `%nin%` <- purrr::compose(`!`, `%in%`)
}
