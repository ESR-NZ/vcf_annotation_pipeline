"""
Author: Miles Benton, Joep de Ligt and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.
Date created: 2020-03-06
Modified: 2020-03-19
Run: snakemake -n -r -j 24 -p --use-conda
Rule diagram: snakemake --rulegraph | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag | dot -Tpng > dag.png
"""

# Adapt paths as appropriate

# Define samples from vcf dir in human_genomics_pipeline using wildcards
SAMPLES, = glob_wildcards("../human_genomics_pipeline/vcf/{sample}.raw.snps.indels.AS.g.vcf")
# Temp dir
TEMPDIR = "/store/lkemp/tmp/"
# links to reference human genome and various annotation databases
GENOME = "/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta"
DBSNP = "/store/lkemp/publicData/dbSNP/ncbi/GRCh37/build151/GATK/All_20180423.vcf.gz"
VEP = "/store/lkemp/publicData/vep/GRCh37/"
DBNSFP = "/store/lkemp/publicData/dbNSFP/GRCh37/dbNSFPv4.0a.hg19.custombuild.gz"
MILLS = "/store/lkemp/publicData/dbSNP/gatkBundle/GRCh37/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
INDEL1000G = "/store/lkemp/publicData/dbSNP/gatkBundle/GRCh37/1000G_phase1.indels.hg19.sites.vcf.gz"
SNP1000G = "/store/lkemp/publicData/dbSNP/gatkBundle/GRCh37/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
OMNI = "/store/lkemp/publicData/genotype/gatkBundle/GRCh37/1000G_omni2.5.hg19.sites.vcf.gz"
HAPMAP = "/store/lkemp/publicData/haplotype/gatkBundle/GRCh37/hapmap_3.3.hg19.sites.vcf.gz"
CADD = "/store/lkemp/publicData/CADD/GRCh37/whole_genome_SNVs.tsv.gz"

rule all:
    input:
        expand("annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf", sample=SAMPLES)

rule gatk4_GenotypeGVCFs:
    input:
        vcf="vcf/{sample}.raw.snps.indels.AS.g.vcf"
    output:
        vcf="vcf/{sample}.genotype.vcf"
    log:
        "logs/gatk_genotype/{sample}.log"
    benchmark:
        "benchmarks/gatk_genotype/{sample}.genotype"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" GenotypeGVCFs \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.vcf} \
            -D {DBSNP} \
            -G StandardAnnotation -G AS_StandardAnnotation
        """

rule gatk4_VariantRecalibrator_indel:
    input:
        vcf="vcf/{sample}.genotype.vcf"
    output:
        report="recalibrated/{sample}.recal.indels",
        tranches="recalibrated/{sample}.tranches.indels",
        rscript="recalibrated/{sample}.plots.indels.R"
    log: 
        "logs/gatk_recal_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_indels/{sample}.recal.indels"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.report} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            --trust-all-polymorphic \
            --max-gaussians 4 \
            -mode INDEL \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {MILLS} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {INDEL1000G} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP}
        """

rule gatk4_VariantRecalibrator_SNP:
    input:
        vcf="vcf/{sample}.genotype.vcf"
    output:
        report="recalibrated/{sample}.recal.snps",
        tranches="recalibrated/{sample}.tranches.snps",
        rscript="recalibrated/{sample}.plots.snps.R"
    log: 
        "logs/gatk_recal_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal_snps/{sample}.snp.recal"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
            -R {GENOME} \
            -V {input.vcf} \
            -O {output.report} \
            --tranches-file {output.tranches} \
            --rscript-file {output.rscript} \
            --trust-all-polymorphic \
            --max-gaussians 6 \
            -mode SNP \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {HAPMAP} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {OMNI} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {SNP1000G} \
            -resource:resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP}
        """

rule gatk4_VQSR_indel:
    input:
        vcf="vcf/{sample}.genotype.vcf",
        recal="recalibrated/{sample}.recal.indels",
        tranches="recalibrated/{sample}.tranches.indels"
    output:
        vcf=temp("recalibrated/{sample}.tmp.vqsr.recal.indels.vcf")
    log:
        "logs/gatk_vqsr_indels/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_indels/{sample}.recal.indels"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk ApplyVQSR -R {GENOME} -V {input.vcf} -O {output.vcf} --recal-file {input.recal} --tranches-file {input.tranches} --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL"

rule gatk4_VQSR_SNP:
    input:
        vcf="recalibrated/{sample}.tmp.vqsr.recal.indels.vcf",
        recal="recalibrated/{sample}.recal.snps",
        tranches="recalibrated/{sample}.tranches.snps"
    output:
        vcf="recalibrated/{sample}.vqsr.recal.vcf"
    log:
        "logs/gatk_vqsr_snps/{sample}.log"
    benchmark:
        "benchmarks/gatk_vqsr_snps/{sample}.recal.snps"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk ApplyVQSR -R {GENOME} -V {input.vcf} -O {output.vcf} --recal-file {input.recal} --tranches-file {input.tranches} --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode SNP"

rule SnpSift:
    input:
        vcf="recalibrated/{sample}.vqsr.recal.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    log: 
        "logs/snpsift/{sample}.log"
    benchmark:
        "benchmarks/snipsift/{sample}.dbnsfp"
    conda:
        "envs/dbnsfp.yaml"
    shell:
        "SnpSift -Xmx16g dbnsfp -v -db {DBNSFP} {input.vcf} > {output.vcf}"

rule VEP:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    log: 
        "logs/vep/{sample}.log"
    benchmark:
        "benchmarks/vep/{sample}.vep"
    conda:
        "envs/vep.yaml"
    threads: 4
    shell:
        "vep -v --assembly GRCh37 --cache --dir {VEP} --fasta {GENOME} -i {input.vcf} -o {output.vcf} --stats_text --everything --vcf --force_overwrite --offline"

rule GENMOD:
    input:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.vcf"
    output:
        vcf="annotated/{sample}.vqsr.recal.dbnsfp.vep.genmod.vcf"
    log: 
        "logs/genmod/{sample}.log"
    benchmark:
        "benchmarks/genmod/{sample}.genmod"
    conda:
        "envs/genmod.yaml"
    threads: 4
    shell:
        "genmod annotate {input.vcf} --regions -c {CADD} -o {output.vcf}"