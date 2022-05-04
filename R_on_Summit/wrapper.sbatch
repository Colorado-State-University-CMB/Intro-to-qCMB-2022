#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:02:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=Rscript
#SBATCH --output=%x.%j.log # gives slurm.ID.log
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

Rscript ./Biological_Databases_mod.Rscript
