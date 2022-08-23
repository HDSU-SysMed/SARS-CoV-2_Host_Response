#!/bin/sh

source ~/.bashrc

conda activate bq

snakemake --cluster "sbatch -t 02:00:00 -n 1 --mem=64gb -J snake" -j 15 --rerun-incomplete 
