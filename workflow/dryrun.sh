#!/bin/bash -x

snakemake \
-n \
-j 32 \
--use-conda \
--conda-frontend mamba \
--latency-wait 20 \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml
