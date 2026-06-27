#!/bin/bash
#SBATCH --job-name=dendro_sil_full
#SBATCH --output=dendro_sil_full_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --cluster=YOUR_CLUSTER
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT

source ~/.bashrc
conda activate thesis

cd /path/to/data/Cell2location/scripts

python 13_2_dendogram_silhouette_non_subsampled.py