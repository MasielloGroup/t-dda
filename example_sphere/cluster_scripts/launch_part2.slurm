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

array=( "785.714_eV" "791.429_eV" "797.143_eV" "802.857_eV" "808.571_eV" "814.286_eV" "820.000_eV" "825.714_eV" "831.429_eV" "837.143_eV" "842.857_eV" "848.571_eV" "854.286_eV" "860.000_eV" "865.714_eV" "871.429_eV" "877.143_eV" "882.857_eV" "888.571_eV" "894.286_eV" "900.000_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



