#!/bin/sh
#SBATCH --job-name='firstjob'
#SBATCH --output=job_%j.out
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=24:00:00
pwd; hostname; date
module load R
Rscript question2code.R