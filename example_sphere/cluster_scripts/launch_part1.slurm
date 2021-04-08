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

array=( "1.536_eV" "1.557_eV" "1.579_eV" "1.600_eV" "1.621_eV" "1.643_eV" "1.664_eV" "1.686_eV" "1.707_eV" "1.729_eV" "1.750_eV" "1.771_eV" "1.793_eV" "1.814_eV" "1.836_eV" "1.857_eV" "1.879_eV" "1.900_eV" "1.921_eV" "1.943_eV" "1.964_eV" "1.986_eV" "2.007_eV" "2.029_eV" "2.050_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



