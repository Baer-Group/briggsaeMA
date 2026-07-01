#!/bin/sh
#SBATCH --job-name=trim # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=40gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=trim_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou


cd /orange/juannanzhou/Rifat_CB_raw/briggsae/Long_read_Fresh

module load hifiadapterfilt
hifiadapterfilt.sh
