#!/bin/bash
# check.sh - check this bundle and the R of the cluster before anything is submitted
#
# Run it on the login node, from the folder of the bundle:
#
#   bash check.sh
#
# It takes seconds, fits nothing and installs nothing. It loads R as settings.sh
# says, then checks that the files of the bundle arrived unchanged, that the
# installed rTFT is the version that wrote the bundle, that cmdstanr and CmdStan
# can be found, and that every model of the bundle can be assembled. Each problem
# it finds is reported with what to do about it, which is most often to run
# install.sh.

set -euo pipefail
cd "$(dirname "$0")"
. ./settings.sh
load_r
if ! Rscript -e 'quit(status = !requireNamespace("rTFT", quietly = TRUE))'; then
  echo "rTFT is not installed in this R. Install it, with what else the bundle needs:"
  echo
  echo "  bash install.sh"
  exit 1
fi
Rscript -e 'rTFT::check_hpc_bundle(".")'
