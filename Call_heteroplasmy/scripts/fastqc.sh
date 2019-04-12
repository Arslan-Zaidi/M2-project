#!/bin/bash

#SBATCH -C new
#SBATCH --cpus-per-task=8
#SBATCH -J fastqc
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=saz5078@psu.edu
#SBATCH -e slurm-%j.err
#SBATCH -D /nfs/brubeck.bx.psu.edu/scratch6/arslan/slrm/M2

source /nfs/brubeck.bx.psu.edu/scratch4/arslan/.bash_profile
source activate cenv

run=${1}
sample=${2}

mkdir -p /nfs/brubeck.bx.psu.edu/scratch3/arslan/TR.rd.trim/qcreps/${run}
fastqc /nfs/brubeck.bx.psu.edu/scratch3/arslan/TR.rd.trim/${run}/${sample}/fil.mapped/nvcReady/srt.${sample}.rd.nvcReady.bam --noextract -o /nfs/brubeck.bx.psu.edu/scratch3/arslan/TR.rd.trim/qcreps/${run} 


