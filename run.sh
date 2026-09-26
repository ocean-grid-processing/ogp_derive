for window in 2005_2024
do
  for level in 0_2000
  do
    sbatch derive.slurm $window $level
  done
done
