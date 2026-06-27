#!/bin/bash
#SBATCH --clusters=wice
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=gsea_lmc
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=gsea_%j.out

source activate thesis
python /path/to/data/14_gsea_lmc.py
