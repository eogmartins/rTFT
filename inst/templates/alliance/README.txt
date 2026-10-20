A bundle of models to fit on a cluster
======================================

Written by {{written_by}} on {{created}}, with rTFT {{rtft_version}}.
{{contents}}

What is in it

  settings.sh   the account, time, memory, cores and modules, in one place
  install.sh    installs rTFT, cmdstanr and CmdStan into the cluster's R
  install.R     the R part of install.sh
  check.sh      checks the bundle and the cluster's R, on the login node
  submit.sh     submits the tasks to Slurm, the cluster's scheduler
  task.sh       what Slurm runs for each task; not run by hand
  status.sh     says what each task has done, and what to do next
  bundle/       the data, the models and the settings of the run, with their
                checksums; nothing in it is to be edited
  stan/         where the compiled Stan programs will be kept
  results/      where the tasks will write what they find
  logs/         where Slurm will write what each task printed

What to do

  1. Copy the folder to the cluster, from your own computer:

       {{copy_command}}

     {{copy_note}}

     Many clusters ask for a second factor of authentication besides your
     password, at each copy and each connection; those of the Digital Research
     Alliance of Canada do. If yours does, it must be set up first, as the
     cluster's documentation explains.

  2. Connect to the cluster and go to the folder:

       ssh {{ssh_target}}
       cd {{remote_folder}}

  3. Install what the bundle needs, once on each cluster (and again for a
     bundle written by another version of rTFT):

       bash install.sh

     The first time, it compiles brms and CmdStan, which can take an hour or
     more. It prints a great deal; its last lines say whether it worked. If
     the login node stops it, run it as a job instead, and wait for the job
     to finish:

       bash install.sh --job

     The job downloads what it installs. On a cluster whose compute nodes
     cannot reach the internet it stops at once and says so; ask the
     cluster's support how R packages are installed there.

  4. Check the bundle and the cluster's R:

       bash check.sh

     check.sh reports anything missing and what to do about it. Run it again
     until it reports nothing missing.

  5. Submit the tasks:

       bash submit.sh

     Each fit is one task of a job array, a set of jobs that Slurm runs side
     by side as cores become free, and each refit for cross-validation is one
     task of a second array, which waits for the fits to end, since each refit
     starts from its model's fit. submit.sh prints the number Slurm gave each
     array. A task can wait in the queue for minutes or for hours, depending
     on how busy the cluster is. The tasks run whether or not you stay
     connected: you can log out, and run status.sh when you come back.

  6. Follow the tasks:

       sq                  the tasks waiting (PD) and running (R)
       bash status.sh      what each task has done, and what to do next

     What each task prints is in logs/, one file per task.

  7. When a task failed, ran out of time or memory, or was stopped by Slurm,
     status.sh says so and gives the command that submits again only the
     tasks not done. If a task ran out of time or memory, raise TIME or
     MEM_PER_CPU in settings.sh first. The command leaves out tasks still
     waiting or running; do not submit those again by hand, or they will run
     twice.

  8. When status.sh says every task is done, copy the results back, in a
     terminal on your own computer, not the one connected to the cluster:

       {{copy_back_command}}

     results/ holds every fit, and can take some gigabytes. If you have moved
     the bundle's folder on your computer since writing it, change the path
     in the command, and in step 9.

  9. Assemble them in R, on your own computer:

       out <- rTFT::collect_hpc("{{local_path}}")
       out$summary

     You can collect before every task is done: models without a result yet
     are named, and step 8 can be run again to bring back the rest.
