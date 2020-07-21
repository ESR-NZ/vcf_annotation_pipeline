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
        expand("annotated/{sample}_filtered_annotated.vcf", sample = SAMPLES)

##### Set up report #####

report: "report/workflow.rst"

##### Load rules #####
if config['DATA'] == "Single" and config['GPU'] == "No":
    include: "rules/gatk_CNNScoreVariants.smk"
    include: "rules/gatk_FilterVariantTranches.smk"
    
if config['DATA'] == "Single" and config['GPU'] == "Yes":
    include: "rules/pbrun_cnnscorevariants.smk"
    include: "rules/gatk_FilterVariantTranches.smk"

if config['DATA'] == "Cohort" and config['GPU'] == "No":
    include: "rules/gatk_VariantRecalibrator_indel.smk"
    include: "rules/gatk_VariantRecalibrator_snp.smk"
    include: "rules/gatk_VQSR_indel.smk"
    include: "rules/gatk_VQSR_snp.smk" # create option for hard filtering if a sample fails vqsr

if config['DATA'] == "Cohort" and config['GPU'] == "Yes":
    include: "rules/pbrun_vqsr.smk" # create option for hard filtering if a sample fails vqsr

include: "rules/SnpSift_dbNSFP.smk"
include: "rules/vep.smk"
include: "rules/genmod_annotate_CADD.smk"