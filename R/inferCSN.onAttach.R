.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("========================================
", pkgname, " version ", version, "

This message can be suppressed by:
  suppressPackageStartupMessages(library(inferCSN))
========================================
")

  packageStartupMessage(inferCSN.logo())
  packageStartupMessage(msg)
}
