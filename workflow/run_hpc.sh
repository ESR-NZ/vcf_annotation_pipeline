#!/bin/bash -x

snakemake \
--cores 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition}"