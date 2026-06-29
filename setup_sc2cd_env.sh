#!/usr/bin/env bash
set -euo pipefail

# Create a conda environment for the SC2CD notebooks.
# Usage:
#   bash setup_sc2cd_env.sh
#   bash setup_sc2cd_env.sh my_env_name

ENV_NAME="${1:-sc2cd_env}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda was not found on PATH. Install Miniconda/Anaconda first." >&2
  exit 1
fi

if conda create --help 2>/dev/null | grep -q -- "--solver"; then
  SOLVER_ARGS=(--solver libmamba)
else
  SOLVER_ARGS=()
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Conda environment '$ENV_NAME' already exists; skipping creation."
else
  echo "Creating conda environment '$ENV_NAME'..."
  conda create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda --override-channels "${SOLVER_ARGS[@]}" \
    python=3.10 \
    ipykernel notebook jupyterlab nbconvert \
    numpy pandas scipy scikit-learn scikit-misc scanpy anndata h5py \
    python-igraph leidenalg louvain \
    matplotlib seaborn statsmodels numba pytables tqdm pot pytorch \
    r-base=4.4 rpy2 \
    r-knitr r-dplyr r-wrmisc r-rcpp r-fields r-gtools \
    r-scatterpie r-ggplot2 r-ggforce r-rcpparmadillo \
    r-ape r-reticulate r-mclust r-rcppml \
    bioconductor-summarizedexperiment bioconductor-singlecellexperiment \
    c-compiler cxx-compiler make
fi

if ! conda list -n "$ENV_NAME" scikit-misc | awk 'NR > 3 {print $1}' | grep -qx "scikit-misc"; then
  echo "Installing missing Scanpy optional dependency scikit-misc..."
  conda install -y -n "$ENV_NAME" -c conda-forge scikit-misc
fi

echo "Creating local writable cache directories..."
mkdir -p .cache/numba .cache/matplotlib .cache/fontconfig

echo "Installing CRAN package NMF inside '$ENV_NAME'..."
conda run -n "$ENV_NAME" Rscript -e '
  env_lib <- file.path(R.home(), "library")
  .libPaths(c(env_lib, setdiff(.libPaths(), env_lib)))
  if (!requireNamespace("NMF", quietly = TRUE)) {
    install.packages("NMF", lib = env_lib, repos = "https://cloud.r-project.org")
  }
'

echo "Registering Jupyter kernel..."
conda run -n "$ENV_NAME" python -m ipykernel install \
  --user \
  --name "$ENV_NAME" \
  --display-name "Python ($ENV_NAME)"

echo "Verifying Python imports..."
conda run -n "$ENV_NAME" env \
  NUMBA_CACHE_DIR="$SCRIPT_DIR/.cache/numba" \
  MPLCONFIGDIR="$SCRIPT_DIR/.cache/matplotlib" \
  XDG_CACHE_HOME="$SCRIPT_DIR/.cache" \
  python -c '
import os
os.environ.setdefault("R_HOME", os.path.join(os.environ["CONDA_PREFIX"], "lib", "R"))
import anndata, h5py, numpy, ot, pandas, scanpy, scipy, sklearn, torch
import skmisc.loess
import rpy2.robjects
from SC2CD_py import SC2CD
import SC2CD_py.CARD_utils
print("Python imports OK")
'

echo "Verifying R packages..."
conda run -n "$ENV_NAME" Rscript -e '
  env_lib <- file.path(R.home(), "library")
  .libPaths(c(env_lib, setdiff(.libPaths(), env_lib)))
  if (is.na(parallel::detectCores())) {
    unlockBinding("detectCores", asNamespace("parallel"))
    assign("detectCores", function(all.tests = FALSE, logical = TRUE) 2L,
           envir = asNamespace("parallel"))
    lockBinding("detectCores", asNamespace("parallel"))
  }
  pkgs <- c(
    "dplyr", "wrMisc", "Rcpp", "fields", "gtools", "scatterpie",
    "ggplot2", "ggforce", "RcppArmadillo", "SummarizedExperiment",
    "SingleCellExperiment", "ape", "NMF", "mclust", "RcppML"
  )
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  cat("R packages OK\n")
'

cat <<EOF

Done.

Next steps:
  conda activate $ENV_NAME
  jupyter notebook

Then select this notebook kernel:
  Python ($ENV_NAME)

EOF

if [[ ! -d Data || ! -f src/CARDfree.cpp || ! -f src/CARDref.cpp ]]; then
  cat <<'EOF'
NOTE:
  The environment is ready, but a full notebook run also needs the project data
  and Rcpp sources:

    Data/
    src/CARDfree.cpp
    src/CARDref.cpp

  Put those in the project root before running the full tutorials.
EOF
fi
