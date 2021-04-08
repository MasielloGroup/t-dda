#!/bin/bash
## Job Name
#SBATCH --job-name=610
## Allocation Definition
#SBATCH --account=chem
#SBATCH --partition=chem
## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (28 is Slurm default)
#SBATCH --ntasks-per-node=28
## Walltime (days-HH:MM:SS)
#SBATCH --time=15:00:00
## Memory per node
#SBATCH --mem=80G
##Output file
#SBTACH --output cfs.out

module load anaconda3_4.3.1

#################################################

array=( "-200" "200" "0" "0" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



