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

array=( "500.000_eV" "505.714_eV" "511.429_eV" "517.143_eV" "522.857_eV" "528.571_eV" "534.286_eV" "540.000_eV" "545.714_eV" "551.429_eV" "557.143_eV" "562.857_eV" "568.571_eV" "574.286_eV" "580.000_eV" "585.714_eV" "591.429_eV" "597.143_eV" "602.857_eV" "608.571_eV" "614.286_eV" "620.000_eV" "625.714_eV" "631.429_eV" "637.143_eV" )

#################################################

for ((i=0;i<${#array[@]};++i)); do
    file=${array[i]}
    cd ${file}
    /gscratch/chem/masiello_group/g-dda/ddscat &> ${file}_DDApump.out &
    cd ..
done; wait



