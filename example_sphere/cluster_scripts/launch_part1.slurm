#!/bin/bash
## Job Name
#SBATCH --job-name=p1
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

array=( "642.857_eV" "648.571_eV" "654.286_eV" "660.000_eV" "665.714_eV" "671.429_eV" "677.143_eV" "682.857_eV" "688.571_eV" "694.286_eV" "700.000_eV" "705.714_eV" "711.429_eV" "717.143_eV" "722.857_eV" "728.571_eV" "734.286_eV" "740.000_eV" "745.714_eV" "751.429_eV" "757.143_eV" "762.857_eV" "768.571_eV" "774.286_eV" "780.000_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



