#!/usr/bin/env bash
set -x

# For the reference: CONDA env vars: https://conda.io/docs/user-guide/tasks/build-packages/environment-variables.html

### Few more R packages that are not on yet conda
#export TAR=/bin/tar  # to avoid "/bin/gtar: not found"
# The "options(unzip = )" hack to address the install_github issue under conda https://github.com/r-lib/devtools/issues/1722
#R -e "library(devtools); options(unzip = '$(which unzip)'); devtools::install_github('Francescojm/CELLector', dependencies=FALSE)"

# Changing permissions to executables
chmod +x ${SRC_DIR}/src/cacao.py
chmod +x ${SRC_DIR}/src/cacao.R
chmod +x ${SRC_DIR}/src/pathogenic_cancer_codons.R
chmod +x ${SRC_DIR}/cacao_wflow.py
# Moving libraries and scripts
mkdir -p ${PREFIX}/lib/R/library/cacao
mkdir -p ${PREFIX}/bin
mv ${SRC_DIR}/src/cacao.py ${SRC_DIR}/cacao_wflow.py ${PREFIX}/bin/  # python scripts
mv ${SRC_DIR}/src/pathogenic_cancer_codons.R ${SRC_DIR}/src/cacao.R ${PREFIX}/bin/  # R scripts
# R modules:
R -e "library(devtools); devtools::install('${SRC_DIR}/src/R/cacao', dependencies=FALSE, args=c('--library=${PREFIX}/lib/R/library'))"
