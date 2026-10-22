#!/bin/bash
# install.sh - install what the bundle needs into the R of the cluster
#
# Run it once on each cluster, and again when the bundle is written by another
# version of rTFT. From the folder of the bundle:
#
#   bash install.sh          on the login node
#   bash install.sh --job    as a Slurm job, when the login node stops a long compile
#
# It loads R as settings.sh says, creates the R library if it does not exist yet, and
# installs what is missing: rTFT at the version that wrote the bundle (with brms and
# the other packages it needs), cmdstanr, and CmdStan, which cmdstanr compiles. Each
# step is skipped when it is already done, so the script can be run again safely.
# Installing brms and CmdStan for the first time compiles a great deal and can take an
# hour or more. When it has finished, run check.sh.

set -euo pipefail
# Slurm runs a copy of the script kept elsewhere, so the install job, known by its
# name, goes back to the folder it was submitted from, which is the bundle's. Run any
# other way, an interactive Slurm session among them, the script goes to the folder
# that holds it.
if [ "${SLURM_JOB_NAME:-}" = "rtft-install" ]; then
  cd "$SLURM_SUBMIT_DIR"
else
  cd "$(dirname "$0")"
fi
. ./settings.sh

if [ "${1:-}" = "--job" ]; then
  mkdir -p logs
  MAIL=()
  if [ -n "$EMAIL" ]; then
    MAIL=(--mail-user="$EMAIL" --mail-type=END,FAIL)
  fi
  sbatch --account="$ACCOUNT" --job-name=rtft-install --time=03:00:00 \
    --cpus-per-task=4 --mem-per-cpu=4000M --output=logs/install-%j.out \
    ${MAIL[@]+"${MAIL[@]}"} "$(basename "$0")"
  echo "Submitted. When the job has finished, read logs/install-<job number>.out,"
  echo "then run:  bash check.sh"
  exit 0
fi

# The compute nodes of some clusters cannot reach the internet, and everything is
# downloaded. The install job finds that out with one request before it installs
# anything, rather than at its first download. (RTFT_CHECK_URL is for rTFT's tests.)
if [ "${SLURM_JOB_NAME:-}" = "rtft-install" ] && command -v curl > /dev/null; then
  if ! curl --silent --head --fail --max-time 30 \
       "${RTFT_CHECK_URL:-https://cloud.r-project.org}" > /dev/null; then
    echo "This job cannot reach the internet: the compute nodes of this cluster have no"
    echo "access to it, so nothing could be downloaded here."
    echo "Run  bash install.sh  on the login node instead, or ask the cluster's support"
    echo "how R packages are installed there."
    exit 1
  fi
fi

# R and the compilers keep their temporary files in TMPDIR, /tmp unless it is set. On
# a login node /tmp is small and shared by everyone connected, and a first install,
# which compiles brms and CmdStan, can fill it. In a Slurm job they go to the node's
# own disk, which Slurm gives each job and clears when it ends; on the login node, to
# the folder tmp/ of the bundle, on the file system that holds the bundle.
if [ -n "${SLURM_TMPDIR:-}" ]; then
  export TMPDIR="$SLURM_TMPDIR"
else
  mkdir -p tmp
  export TMPDIR="$PWD/tmp"
fi

load_r
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-2}"
Rscript install.R
echo
echo "Done. Now run:  bash check.sh"
