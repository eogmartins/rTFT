# settings.sh - the one file where the settings of this bundle are written
#
# Written by rTFT::write_hpc_bundle() on {{created}}, for the clusters of the
# Digital Research Alliance of Canada. Every script of the bundle reads this file,
# so a value changed here changes for all of them. Make sure to change it here and
# nowhere else. Change the time or the memory here when a task ran out of either,
# and submit again.

# The allocation the jobs are charged to, as `sbatch --account` takes it.
ACCOUNT="{{account}}"

# Where Slurm sends an email when a task fails; empty for no email.
EMAIL="{{email}}"

# The time each task may run (hours:minutes:seconds, or days-hours:minutes:seconds),
# the memory per core, and the cores per task, one for each chain. The time is the
# task's wall time, the time it may run by the clock; each cluster sets a longest
# wall time a job may ask for, and Slurm refuses a task that asks for more.
TIME="{{time}}"
MEM_PER_CPU="{{mem_per_cpu}}"
CPUS_PER_TASK="{{cpus_per_task}}"

# The modules that give R on this cluster. `module avail r` lists the versions of R
# that the cluster offers.
MODULES="{{modules}}"

# The R library that holds rTFT, brms and cmdstanr; empty for R's own personal
# library.
R_LIBRARY="{{r_library}}"

# Load R as the settings above say, from a clean set of modules. The clusters of the
# Alliance keep some modules loaded through a purge, and say so at length each time;
# that notice is not shown. Errors in loading the modules below still are.
load_r() {
  module purge > /dev/null 2>&1
  module load $MODULES
  if [ -n "$R_LIBRARY" ]; then
    export R_LIBS="$R_LIBRARY"
  fi
}
