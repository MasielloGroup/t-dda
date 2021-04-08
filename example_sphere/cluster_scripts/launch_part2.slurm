#!/bin/bash
## Job Name
#SBATCH --job-name=p2
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

array=( "2.071_eV" "2.093_eV" "2.114_eV" "2.136_eV" "2.157_eV" "2.179_eV" "2.200_eV" "2.221_eV" "2.243_eV" "2.264_eV" "2.286_eV" "2.307_eV" "2.329_eV" "2.350_eV" "2.371_eV" "2.393_eV" "2.414_eV" "2.436_eV" "2.457_eV" "2.479_eV" "2.500_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



