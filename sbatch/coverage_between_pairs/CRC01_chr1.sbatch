#!/bin/bash
#SBATCH --output=/cs/usr/drordod/Desktop/project/proj_scwgbs/output/coverage_between_pairs_CRC01_chr1.out
#SBATCH --error=/cs/usr/drordod/Desktop/project/proj_scwgbs/output/coverage_between_pairs_CRC01_chr1.err
#SBATCH --time=72:00:00
#SBATCH --mem=4G

RUNPATH=~/Desktop/project/proj_scwgbs
cd $RUNPATH/format_files
python3 coverage_between_pairs.py --cpg_format_files /cs/usr/drordod/Desktop/scTrio-seq-reanalysis/liordror/cpg_format/filtered/CRC01/all_cpg_ratios_CRC01_chr1.dummy.pkl.zip --output_folder /cs/usr/drordod/Desktop/scTrio-seq-reanalysis/liordror/stats/coverage_histograms/
