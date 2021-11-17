#!/bin/bash
## Job Name
#SBATCH --job-name=610
## Allocation Definition
#SBATCH --account=chem-ckpt
#SBATCH --partition=ckpt
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

/gscratch/chem/masiello_group/g-dda/source_code/ddscat
wait
python sphere_qabs.py
wait
/gscratch/chem/masiello_group/t-dda/source_code/Lattice_Diffusion /gscratch/chem/masiello_group/myGreen.num_300 var.par tdda_input temp.out



