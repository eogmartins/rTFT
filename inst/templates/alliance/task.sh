#!/bin/bash
# task.sh - one task of a job array: one fit, or one refit for cross-validation
#
# submit.sh hands this script to Slurm, which runs it once for each task of the
# array, on a compute node, with the number of the task in SLURM_ARRAY_TASK_ID. It
# loads R as settings.sh says and calls rTFT::run_hpc_task(), which fits the model of
# that task, or refits it without one fold, and writes the result into results/.
# It is not run by hand.

set -euo pipefail
# Slurm runs a copy of the script kept elsewhere, so the task goes back to the folder
# it was submitted from, which is the bundle's.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
. ./settings.sh

STAGE="${1:-}"
if [ "$STAGE" != "fit" ] && [ "$STAGE" != "kfold" ]; then
  echo "task.sh takes the stage, fit or kfold."
  exit 1
fi
echo "Task ${SLURM_ARRAY_TASK_ID:-} of the $STAGE stage, on $(hostname), $(date)"
load_r
# R's temporary files, the draws cmdstanr writes for each chain among them, go to the
# node's own disk, which Slurm gives each job and clears when it ends, and not to
# /tmp, which every job on the node shares; the fit holds its draws before it is
# saved, so nothing is lost when the disk is cleared.
export TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
Rscript -e "rTFT::run_hpc_task('.', '$STAGE')"
