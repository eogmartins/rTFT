#!/bin/bash
# status.sh - what each task of this bundle has done, and what to do next
#
# Run it on the login node, from the folder of the bundle, at any time after
# submit.sh:
#
#   bash status.sh
#
# It loads R as settings.sh says and calls rTFT::hpc_status(), which reads what the
# tasks have written into results/ and logs/. It says how many tasks of each stage
# are done, which are not and why, and gives the command that submits again the
# tasks not done, leaving out those still waiting or running. It changes nothing.

set -euo pipefail
cd "$(dirname "$0")"
. ./settings.sh
load_r
Rscript -e 'invisible(rTFT::hpc_status("."))'
