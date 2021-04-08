#!/bin/bash
## Job Name
#SBATCH --job-name=1t25
## Allocation Definition
#SBATCH --account=chem-ckpt
#SBATCH --partition=ckpt
## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (28 is Slurm default)
#SBATCH --ntasks-per-node=28
## Walltime (days-HH:MM:SS)
#SBATCH --time=48:00:00
## Memory per node
#SBATCH --mem=100G
##Output file
#SBTACH --output cfs.out

for ((i=1;i<=25;i++));
    do 
    cd W$i
    ln -s /gscratch/chem/masiello_group/srcPW_background/ddscat 
    ./ddscat &> W$i.out &
    cd ..
    done
wait
       
