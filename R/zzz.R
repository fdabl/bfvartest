.onLoad <- function(libname, pkgname) {
  if (length(stanmodels) != 0) {
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) loadModule(m, what = TRUE)
  } else {
    message("No stan programs to compile were found.")
  }
}
