#!/bin/bash
#SBATCH --clusters=wice
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=hierarchical_clustering
#SBATCH --time=2:00:00
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/path/to/data/Cell2location/spatial_mapping/exploration_mean/hierarchical_clustering_%j.log
#SBATCH --error=/path/to/data/Cell2location/spatial_mapping/exploration_mean/hierarchical_clustering_%j.err

source ~/.bashrc
conda activate thesis

python /path/to/data/Cell2location/scripts/13_run_hierarchical.py
