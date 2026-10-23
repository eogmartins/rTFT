#!/bin/bash
# submit.sh - submit the tasks of this bundle to Slurm
#
# Run it on the login node, from the folder of the bundle, once check.sh has passed:
#
#   bash submit.sh                          every task: the fits, then the refits
#   bash submit.sh --fit 2,5                only fits 2 and 5
#   bash submit.sh --fit 2 --kfold 3-4      fit 2, and refits 3 and 4 after it
#
# Each stage is one job array, with one task for each row of its table in inputs/,
# and each task is given the time, memory and cores of settings.sh. A refit for
# cross-validation starts from its model's fit, so the refits are queued to start
# only when every fit submitted with them has ended, whether it ended well or not; a
# refit whose fit failed stops at once and says so. Submitting again a task that has
# ended well costs seconds, as it is not run again; a task still waiting or running
# would run twice. status.sh gives the tasks that need submitting, and leaves out
# those still waiting or running.

set -euo pipefail
cd "$(dirname "$0")"
. ./settings.sh

usage() {
  echo "Usage: bash submit.sh [--fit TASKS] [--kfold TASKS]"
  echo "TASKS are numbers from 1, as Slurm's --array takes them, such as 3 or 1-4,7."
  exit 1
}

# the number of tasks of a stage: the rows of its table, without the header
count() { echo $(( $(wc -l < "$1") - 1 )); }
N_FIT=$(count inputs/tasks-fit.csv)
N_KFOLD=0
if [ -f inputs/tasks-kfold.csv ]; then
  N_KFOLD=$(count inputs/tasks-kfold.csv)
fi

FIT=""
KFOLD=""
if [ $# -eq 0 ]; then
  FIT="1-$N_FIT"
  if [ "$N_KFOLD" -gt 0 ]; then
    KFOLD="1-$N_KFOLD"
  fi
fi
while [ $# -gt 0 ]; do
  if [ $# -lt 2 ]; then usage; fi
  case "$1" in
    --fit)   FIT="$2" ;;
    --kfold) KFOLD="$2" ;;
    *)       usage ;;
  esac
  if ! [[ "$2" =~ ^[1-9][0-9]*(-[1-9][0-9]*)?(,[1-9][0-9]*(-[1-9][0-9]*)?)*$ ]]; then
    usage
  fi
  shift 2
done
if [ -n "$KFOLD" ] && [ "$N_KFOLD" -eq 0 ]; then
  echo "This bundle has no refits for cross-validation:"
  echo "it was written with kfold = FALSE."
  exit 1
fi
# A task beyond the table would be submitted, and stop at once without a record.
largest() { echo "$1" | tr ',-' '\n\n' | sort -n | tail -1; }
if [ -n "$FIT" ] && [ "$(largest "$FIT")" -gt "$N_FIT" ]; then
  echo "This bundle has $N_FIT fits; fit $(largest "$FIT") is not one of them."
  exit 1
fi
if [ -n "$KFOLD" ] && [ "$(largest "$KFOLD")" -gt "$N_KFOLD" ]; then
  echo "This bundle has $N_KFOLD refits; refit $(largest "$KFOLD") is not one of them."
  exit 1
fi

MAIL=()
if [ -n "$EMAIL" ]; then
  MAIL=(--mail-user="$EMAIL" --mail-type=FAIL)
fi
# Record that the tasks $2 of the stage $1 were submitted as the job $3, with the
# time by the clock, in results/submitted.csv. hpc_status() reads it to tell a task
# waiting in the queue from one that stopped, and does not give it to submit again.
record() {
  if [ ! -f results/submitted.csv ]; then
    echo '"time","stage","tasks","job"' > results/submitted.csv
  fi
  echo "$(date +%s),\"$1\",\"$2\",\"$3\"" >> results/submitted.csv
}

# Submit the tasks $2 of the stage $1 as one job array, with any further options of
# sbatch, and print the number Slurm gives the job.
submit() {
  local stage="$1" tasks="$2"
  shift 2
  sbatch --parsable --account="$ACCOUNT" --job-name="rtft-$stage" --array="$tasks" \
    --time="$TIME" --mem-per-cpu="$MEM_PER_CPU" --cpus-per-task="$CPUS_PER_TASK" \
    --output="logs/$stage-%A_%a.out" ${MAIL[@]+"${MAIL[@]}"} "$@" task.sh "$stage" |
    cut -d ';' -f 1
}

mkdir -p logs results stan
AFTER=()
if [ -n "$FIT" ]; then
  FIT_JOB=$(submit fit "$FIT")
  record fit "$FIT" "$FIT_JOB"
  echo "Submitted the fits $FIT of $N_FIT as job $FIT_JOB."
  AFTER=(--dependency="afterany:$FIT_JOB")
fi
if [ -n "$KFOLD" ]; then
  KFOLD_JOB=$(submit kfold "$KFOLD" ${AFTER[@]+"${AFTER[@]}"})
  record kfold "$KFOLD" "$KFOLD_JOB"
  if [ -n "$FIT" ]; then
    echo "Submitted the refits $KFOLD of $N_KFOLD as job $KFOLD_JOB,"
    echo "to start when the fits have ended."
  else
    echo "Submitted the refits $KFOLD of $N_KFOLD as job $KFOLD_JOB."
  fi
fi
echo
echo "Follow them with:  sq               the tasks waiting (PD) and running (R)"
echo "                   bash status.sh   what each task has done, and what to do next"
