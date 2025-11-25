#!/bin/bash
#SBATCH -A Research_Project-T112069 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D /lustre/projects/Research_Project-T112069/Meth/QTL/Scripts # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=j.d.harvey@exeter.ac.uk # email address
#SBATCH --error=/lustre/projects/Research_Project-T112069/Meth/QTL/Scripts/errorFile.txt

cd /lustre/projects/Research_Project-T112069/Meth/QTL/Data/
source $1


module load R/4.2.1-foss-2022a
module load PLINK/1.07-foss-2016a 

OutDir=$(dirname "$OutPrefix")
mkdir -p $OutDir

IFS=',' read -r -a array <<< "$chr"

if [ $chr == all ]
then
  array=($(seq 1 1 22))  ## seq FIRST STEP LAST
fi

for i in "${array[@]}"
do
  if [ ! -f ${OutPrefix}.chr${i}.raw ]
  then
    echo "Running chr ${i}..."
    plink --bfile ${GenotypeBinaryPrefix} --recodeA --chr $i --out ${OutPrefix}.chr${i}
    echo "#########################################################################################"
  fi
done

#if [ ! -f ${OutPrefix}.eigenvec ]
#then
#  plink --bfile ${GenotypeBinaryPrefix} --indep-pairwise 50 5 0.2 --out pruned
#  plink --bfile ${GenotypeBinaryPrefix} --extract pruned.prune.in --pca --out ${OutPrefix}
#fi 


Rscript /lustre/projects/Research_Project-T112069/Meth/QTL/Scripts/QTL-Analysis/PrepareData.R ${GenotypeBinaryPrefix}.fam ${OutPrefix}.eigenvec $ExpressionFile $PhenotypeFile $GeneLocationFile "$FactCovar" "$NumCovar" "$OutPrefix" $chr

