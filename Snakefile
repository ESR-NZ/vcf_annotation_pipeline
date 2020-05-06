"""
Author: Miles Benton, Joep de Ligt and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.
Date created: 2020-03-06
Modified: 2020-04-29
Dry run: snakemake -n -j 24 --use-conda --use-singularity --singularity-args '-B /dir/to/databases/' --configfile config.yaml
Full run: snakemake -j 24 --use-conda --use-singularity --singularity-args '-B /dir/to/databases/' --configfile config.yaml
Report: snakemake --report report.html --configfile config.yaml --report-stylesheet stylesheet.css
Rule diagram: snakemake --rulegraph --configfile config.yaml | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag --configfile config.yaml | dot -Tpng > dag.png
"""

##### Set up #####

# Define samples from vcf dir using wildcards
SAMPLES, = glob_wildcards("../vcf/{sample}_raw_snps_indels_AS_g.vcf")

##### Target rules #####

rule all:
    input:
        expand("annotated/{sample}_filtered_dbnsfp_vep.vcf_summary.txt", sample = SAMPLES),
        expand("annotated/{sample}_filtered_dbnsfp_vep.vcf_warnings.txt", sample = SAMPLES),
        expand("annotated/{sample}_filtered_annotated.vcf", sample = SAMPLES)

##### Set up report #####

report: "report/workflow.rst"

##### Load rules #####
if config['DATA'] == "Single":
    include: "rules/gatk_cnn_score_variants.smk"
    include: "rules/gatk_filter_variant_tranches.smk"
elif config['DATA'] == "Cohort":
    include: "rules/gatk_variant_recalibrator_indel.smk"
    include: "rules/gatk_variant_recalibrator_snp.smk"
    include: "rules/gatk_vqsr_indel.smk"
    include: "rules/gatk_vqsr_snp.smk"
else:
    print: ("ERROR: Please check the values you provided in the configuration file")

include: "rules/snpsift_dbnsfp.smk"
include: "rules/vep.smk"
include: "rules/genmod_cadd.smk"