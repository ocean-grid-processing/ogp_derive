#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=ohc_derive
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=256GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
# MLD: a single submission, so no --level/--levels (the identity level); the field as published, with
# its ensemble spread. as_published verifies the footprint our publish step enforced.
#   sbatch derive_mld.slurm
tag=260921-OP20260507
product_name=LocalGP
author=Giglio_etal2026
citation="Giglio, D. et al. (2026): LocalGP mixed layer depth. DOI: <fill in on deposit>."
results=/scratch/alpine/wimi7695/mld_post/results-260921
python run.py \
    ${results}/MLD_${tag}*.nc \
    --quantities field \
    --bathy /scratch/alpine/wimi7695/ohc_prod/params/etopo60.cdf \
    --mask as_published \
    --tag  ${tag} \
    --provenance-link https://github.com/ocean-grid-processing/provenance/blob/main/README.md \
    --code-version https://github.com/ocean-grid-processing/ohc_derive/releases/tag/1.1.0 \
    --product-name ${product_name} --author ${author} --citation "${citation}" \
    --out ${results}
