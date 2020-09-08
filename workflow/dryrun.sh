#!/bin/bash -x

snakemake \
-n \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml