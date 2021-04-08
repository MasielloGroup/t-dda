#!/bin/bash
## Job Name
#SBATCH --job-name=p0
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

array=( "1.000_eV" "1.021_eV" "1.043_eV" "1.064_eV" "1.086_eV" "1.107_eV" "1.129_eV" "1.150_eV" "1.171_eV" "1.193_eV" "1.214_eV" "1.236_eV" "1.257_eV" "1.279_eV" "1.300_eV" "1.321_eV" "1.343_eV" "1.364_eV" "1.386_eV" "1.407_eV" "1.429_eV" "1.450_eV" "1.471_eV" "1.493_eV" "1.514_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



