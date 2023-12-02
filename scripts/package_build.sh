#!/bin/bash

# If add new C/C++ function(s)
Rscript -e 'Rcpp::compileAttributes()'

# Build new document
Rscript -e 'devtools::document()'

# Fast check test and example(s)
Rscript -e 'devtools::test()'
Rscript -e 'devtools::run_examples()'

# Check R package
Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))'

# check on other distributions
# _win devel
Rscript -e 'devtools::check_win_devel()'

# Verify you're ready for release, and release
Rscript -e 'devtools::release()'
