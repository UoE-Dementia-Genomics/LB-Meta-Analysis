#!/bin/bash
#SBATCH -A Research_Project-T112069 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D /lustre/projects/Research_Project-T112069/Meth/QTL/Data # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=2:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=j.d.harvey@exeter.ac.uk # email address

module load R/4.2.1-foss-2022a

Rscript $1/Main.R $2 $SLURM_ARRAY_TASK_ID $3 $4 $5 $6 

