.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("========================================
", pkgname, " version ", version, "

This message can be suppressed by:
  suppressPackageStartupMessages(library(inferCSN))

To restore the default display use options(\"restore_inferCSN_show\" = TRUE)
========================================
")

  packageStartupMessage(inferCSN.logo())
  packageStartupMessage(msg)
}
