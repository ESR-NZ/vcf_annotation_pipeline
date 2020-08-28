#!/bin/bash -x

snakemake \
-j 32 \
--use-conda \
--conda-frontend mamba \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml