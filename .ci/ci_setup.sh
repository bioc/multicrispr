#!/bin/bash

set -eo pipefail

# install BiocManager
Rscript -e 'install.packages("BiocManager", repos="http://cran.r-project.org")'
# update BiocManager to devel
Rscript -e 'BiocManager::install(version = "devel", update = TRUE, ask = FALSE)'
# reinstall all packages to avoid "Error: package X was installed before R 4.0.0: please re-install it"
Rscript -e 'BiocManager::install(installed.packages()[,1])'
# install r dependencies
# https://www.cyberciti.biz/faq/unix-howto-read-line-by-line-from-file/
while IFS= read -r package;
do
  RDscript -e "BiocManager::install('"$package"', update = TRUE, ask = FALSE)";
done < ".ci/r-requirements.txt"
