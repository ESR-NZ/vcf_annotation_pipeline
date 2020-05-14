#!/bin/bash -x

snakemake \
-j 24 \
--use-conda \
--use-singularity \
--singularity-args '-B /dir/to/databases/' \
--configfile config.yaml