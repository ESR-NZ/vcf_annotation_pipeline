"""
Author: Miles Benton, Joep de Ligt and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.
Date created: 2020-03-06
Modified: 2020-03-19
Run: snakemake -n -r -j 24 -p --use-conda --use-singularity --configfile config_GRCh37.yaml
Rule diagram: snakemake --rulegraph --configfile config_GRCh37.yaml | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag --configfile config_GRCh37.yaml | dot -Tpng > dag.png
"""

##### Set up #####

# Define samples from vcf dir in human_genomics_pipeline using wildcards
SAMPLES,=glob_wildcards("../vcf/{sample}.raw.snps.indels.AS.g.vcf")

if config['BUILD'] == "GRCh37":
    MODELRULE = "rules/model_37.smk"
elif config['BUILD'] == "GRCh38":
    MODELRULE = "rules/model_38.smk"
else: print("ERROR: Please choose either the GRCh37 or GRCh38 build of the reference human genome")

##### Target rules #####

rule all:
    input:
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf", sample = SAMPLES)

##### Set up report #####

report: config["REPORTWORKFLOW"]

##### Load rules #####

include: "rules/genotype.smk"
include: MODELRULE
include: "rules/recalibrate.smk"
include: "rules/annotate.smk"
