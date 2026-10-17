A bundle of models to fit on a cluster
======================================

Written by rTFT::write_hpc_bundle() on {{created}}, with rTFT {{rtft_version}}.
{{contents}}

What is in it

  settings.sh   the account, time, memory, cores and modules, in one place
  install.sh    installs rTFT, cmdstanr and CmdStan into the cluster's R
  install.R     the R part of install.sh
  check.sh      checks the bundle and the cluster's R, on the login node
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
