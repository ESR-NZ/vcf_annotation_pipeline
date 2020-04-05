"""
Author: Miles Benton, Joep de Ligt and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.
Date created: 2020-03-06
Modified: 2020-03-19
Run: snakemake -n -r -j 24 -p --use-conda --use-singularity
Rule diagram: snakemake --rulegraph | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag | dot -Tpng > dag.png
"""

##### load config and set global variables #####

configfile: "config.yaml"

PUBLICDIR = config["PUBLICDIR"]

SAMPLEDIR = config["SAMPLEDIR"]

# Define samples from vcf dir in human_genomics_pipeline using wildcards
SAMPLES, = glob_wildcards("../human_genomics_pipeline/vcf/{sample}.raw.snps.indels.AS.g.vcf")

# links to reference human genome and various annotation databases
GENOME = expand("{publicdir}ucsc.hg19.fasta", publicdir=PUBLICDIR)
DBSNP = expand("{publicdir}All_20180423.vcf.gz", publicdir=PUBLICDIR)
VEP = expand("{publicdir}GRCh37/", publicdir=PUBLICDIR)
DBNSFP = expand("{publicdir}dbNSFPv4.0a.hg19.custombuild.gz", publicdir=PUBLICDIR)
MILLS = expand("{publicdir}Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz", publicdir=PUBLICDIR)
INDEL1000G = expand("{publicdir}1000G_phase1.indels.hg19.sites.vcf.gz", publicdir=PUBLICDIR)
SNP1000G = expand("{publicdir}1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz", publicdir=PUBLICDIR)
OMNI = expand("{publicdir}1000G_omni2.5.hg19.sites.vcf.gz", publicdir=PUBLICDIR)
HAPMAP = expand("{publicdir}hapmap_3.3.hg19.sites.vcf.gz", publicdir=PUBLICDIR)
CADD = "../../vcf_annotation_pipeline/CADD/whole_genome_SNVs.tsv.gz"

##### target rules #####

rule all:
    input:
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf", sample=SAMPLES)

##### load rules #####

include: "rules/genotype.smk"
include: "rules/recalibrate.smk"
include: "rules/annotate.smk"