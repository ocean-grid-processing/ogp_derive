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
# OHC: one synthetic level ($2) from the LocalGP plan, one baseline window ($1); run.sh loops both.
#   sbatch derive_ohc.slurm 2005_2024 0_2000
tag=260908-OP20260507
product_name=LocalGP
author=Giglio_etal2026
citation="Giglio, D. et al. (2026): LocalGP ocean heat content — mask/coverage auxiliaries. DOI: <fill in on deposit>."
results=/scratch/alpine/wimi7695/ohc_prod/results_260908
python run.py \
    ${results}/OHC_${tag}*.nc \
    --levels levels/localgp.toml --level $2 \
    --time-window $1 \
    --quantities ohca,ohu,ohca_trend,ohu_trend,map \
    --mask contiguous_from_top \
    --bathy /scratch/alpine/wimi7695/ohc_prod/params/etopo60.cdf \
    --tag  ${tag} \
    --provenance-link https://github.com/ocean-grid-processing/provenance/blob/main/README.md \
    --code-version https://github.com/ocean-grid-processing/ohc_derive/releases/tag/1.1.0 \
    --product-name ${product_name} --author ${author} --citation "${citation}" \
    --out ${results}
