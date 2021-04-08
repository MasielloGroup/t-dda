for file in *.slurm;
do
  sbatch $file
done
