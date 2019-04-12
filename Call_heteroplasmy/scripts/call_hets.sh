#!/bin/bash

#SBATCH -C new
#SBATCH --cpus-per-task=8
#SBATCH -J hetcall
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=saz5078@psu.edu
#SBATCH -e slurm-%j.err
#SBATCH -D /nfs/brubeck.bx.psu.edu/scratch6/arslan/slrm/M2

source /nfs/brubeck.bx.psu.edu/scratch4/arslan/.bash_profile
source activate cenv

wdir=${1}
cfile=${2}
ofile=${3}
maft=${4}
sdepth=${5}

cd ${wdir}
call_heteroplasmy.py -c ${cfile} -o ${ofile} -d ${sdepth} -f ${maft}
