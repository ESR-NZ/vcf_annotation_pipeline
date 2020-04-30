"""
Author: Miles Benton, Joep de Ligt and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.
Date created: 2020-03-06
Modified: 2020-04-29
Dry run: snakemake -n -j 24 --use-conda --use-singularity --singularity-args '-B /your/dir/' --configfile your_config.yaml
Full run: snakemake -j 24 --use-conda --use-singularity --singularity-args '-B /your/dir/' --configfile your_config.yaml
Rule diagram: snakemake --rulegraph --use-singularity --configfile your_config.yaml | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag --configfile your_config.yaml | dot -Tpng > dag.png
"""

##### Set up #####

# Define samples from vcf dir in human_genomics_pipeline using wildcards
SAMPLES,=glob_wildcards("../vcf/{sample}.raw.snps.indels.AS.g.vcf")

# Define which variation of the modelling rules should used based on the build of reference human genome
# Also define which overall workflow description should be used in the final report
if config['BUILD'] == "GRCh37":
    MODELRULE = "rules/model_37.smk"
    REPORTWORKFLOW = "report/workflow_37.rst"
elif config['BUILD'] == "GRCh38":
    MODELRULE = "rules/model_38.smk"
    REPORTWORKFLOW = "report/workflow_38.rst"
else: print("ERROR: Please choose either the GRCh37 or GRCh38 build of the reference human genome")

##### Target rules #####

rule all:
    input:
        expand("recalibrated/{sample}.plots.indels.R.pdf", sample = SAMPLES),
        expand("recalibrated/{sample}.plots.snps.R.pdf", sample = SAMPLES),
        expand("recalibrated/{sample}.tranches.snps.pdf", sample = SAMPLES),
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf_summary.txt", sample = SAMPLES),
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf_warnings.txt", sample = SAMPLES),
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf", sample = SAMPLES)

##### Set up report #####

report: REPORTWORKFLOW

##### Load rules #####

include: "rules/genotype.smk"
include: MODELRULE
include: "rules/recalibrate.smk"
include: "rules/annotate.smk"
