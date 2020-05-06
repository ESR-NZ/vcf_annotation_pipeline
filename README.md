# vcf_annotation_pipeline

A simple Snakemake workflow to annotate variant call format (VCF) files using GATK4, SnpSift, VEP and genmod. Designed to be used after human_genomics_pipeline.

- [vcf_annotation_pipeline](#vcfannotationpipeline)
  - [Workflow diagram](#workflow-diagram)
  - [How to run vcf_annotation_pipeline](#how-to-run-vcfannotationpipeline)
    - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
    - [2. Take the pipeline to the data on your local machine](#2-take-the-pipeline-to-the-data-on-your-local-machine)
    - [3. Create a local copy of vcf annotation databases (either GRCh37 or GRCh38)](#3-create-a-local-copy-of-vcf-annotation-databases-either-grch37-or-grch38)
      - [GRCh37](#grch37)
      - [GRCh38](#grch38)
    - [4. Choose and modify an appropriate configuration file](#4-choose-and-modify-an-appropriate-configuration-file)
    - [5. Create and activate a conda environment with python and snakemake and installed](#5-create-and-activate-a-conda-environment-with-python-and-snakemake-and-installed)
    - [6. Run the pipeline](#6-run-the-pipeline)
    - [7. Evaluate the pipeline run](#7-evaluate-the-pipeline-run)
    - [8. Commit and push to your forked version of the repo](#8-commit-and-push-to-your-forked-version-of-the-repo)
    - [9. Create a pull request with the upstream repo to merge any useful changes (optional)](#9-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-optional)
  - [Useful reading](#useful-reading)

## Workflow diagram

<img src="rulegraph.png" class="center">

## How to run vcf_annotation_pipeline

- **Prerequisite software:** [R](https://www.r-project.org/), [Git](https://git-scm.com/), [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/), [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/), [samtools](http://www.htslib.org/)
- **Prerequisite data:** Reference human genome and dbSNP database (see [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) for more information on downloading this data)
- **OS:** Validated on Ubuntu 16.04

### 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help forking a github repository

### 2. Take the pipeline to the data on your local machine

Clone the forked [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) repo into the same directory as your vcf data to be processed. Required folder structure and file naming convention:

```bash

.
|___vcf/
|     |___sample1_raw_snps_indels_AS_g.vcf
|     |___sample2_raw_snps_indels_AS_g.vcf
|     |___ ...
|
|___bams/
|     |___sample1_bwa_recal.bam
|     |___sample2_bwa_recal.bam
|     |___ ...
|
|___vcf_annotation_pipeline/

```

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help cloning a forked github repository.

### 3. Create a local copy of vcf annotation databases (either GRCh37 or GRCh38)

#### GRCh37

Download [Ensembl-VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) database using a [conda version of Ensembl-VEP](https://anaconda.org/bioconda/ensembl-vep)

```bash
conda create -n download_data_env python=3.7
conda activate download_data_env
conda install -c bioconda ensembl-vep=99.2
vep_install -a cf -s homo_sapiens -y GRCh37 -c /output/file/path/GRCh37 --CONVERT
conda deactivate
```

Download other databases from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
```

Convert these files to bgzip format in order to create their associated index files

```bash
gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
bgzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

gunzip 1000G_phase1.indels.hg19.sites.vcf.gz
bgzip 1000G_phase1.indels.hg19.sites.vcf
tabix 1000G_phase1.indels.hg19.sites.vcf.gz

gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
bgzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf
tabix 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

gunzip 1000G_omni2.5.hg19.sites.vcf.gz
bgzip 1000G_omni2.5.hg19.sites.vcf
tabix 1000G_omni2.5.hg19.sites.vcf.gz

gunzip hapmap_3.3.hg19.sites.vcf.gz
bgzip hapmap_3.3.hg19.sites.vcf
tabix hapmap_3.3.hg19.sites.vcf.gz
```

Download the [CADD database](https://cadd.gs.washington.edu/download) and it's associated index file.

```bash
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
```

Create a custom [dbNSFP database](https://sites.google.com/site/jpopgen/dbNSFP) build by following [this documentation](https://github.com/GenomicsAotearoa/dbNSFP_build)

#### GRCh38

Download [Ensembl-VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) database using a [conda version of Ensembl-VEP](https://anaconda.org/bioconda/ensembl-vep)

```bash
conda create -n download_data_env python=3.7
conda activate download_data_env
conda install -c bioconda ensembl-vep=99.2
vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/file/path/GRCh38 --CONVERT
conda deactivate
```

Download other databases from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz.tbi
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz.tbi
```

Download the [CADD database](https://cadd.gs.washington.edu/download) and it's associated index file.

```bash
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

Create a custom [dbNSFP database](https://sites.google.com/site/jpopgen/dbNSFP) build by following [this documentation](https://github.com/GenomicsAotearoa/dbNSFP_build)

### 4. Choose and modify an appropriate configuration file

Choose the appropriate config file:

- Use 'GRCh37_config_template.yaml' to run the pipeline against the GRCh37 reference genome
- Use 'GRCh38_config_template.yaml' to run the pipeline against the GRCh38 reference genome

Ensure this choice is defined in your configuration file. For example:

```yaml
# Specify the build of reference genome used (either 'GRCh37' or 'GRCh38')
BUILD: "GRCh38"
```

Specify whether the data is to be analysed on it's own ('Single') or as a part of a cohort ('Cohort'). For example:

```yaml
# Specify the type of input data (either 'Single' or 'Cohort')
DATA: "Single"
```

Set the the working directories in the config file to the reference human genome file, dbSNP database file and the various vcf annotation database files (GRCh37 or GRCh38). For example:

*Note. there is no 1000G indel file available for the GRCh38 build of the reference human genome*

```yaml
# File directories to reference genome, dbSNP database and various vcf annotation databases
REFGENOME: "/home/lkemp/publicData/referenceGenome/Homo_sapiens_assembly38.fasta.gz"
dbSNP: "/home/lkemp/publicData/dbSNP/All_20180418.vcf.gz"
VEP: "/home/lkemp/publicData/VEP/GRCh38/"
dbNSFP: "/home/lkemp/publicData/dbNSFP/dbNSFPv4.0a_custombuild.gz"
MILLS: "/home/lkemp/publicData/Mills/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
SNP1000G: "/home/lkemp/publicData/1000G/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
OMNI: "/home/lkemp/publicData/1000G/1000G_omni2.5.hg38.vcf.gz"
HAPMAP: "/home/lkemp/publicData/hapmap/hapmap_3.3.hg38.vcf.gz"
CADD: "/home/lkemp/publicData/CADD/whole_genome_SNVs.tsv.gz"
```

If analysing WES data, pass a design file (.bed) indicating the genomic regions that were sequenced (see [here](https://leahkemp.github.io/documentation/human_genomic_pipelines/design_files.html) for more information on accessing design files). Also set the level of padding. For example:

*Note: If not analysing WES data, leave these fields blank*

```yaml
# Whole exome sequencing (WES) specific settings (leave blank if analysing other data such as WGS)
WES:
  # Genomic intervals over which to operate
  INTERVALS: "-L /home/lkemp/publicData/sure_select_human_all_exon_V7/S31285117_AllTracks.bed"
  # Amount of padding (in bp) to add to each interval you are including
  PADDING: "-ip 100"
```

If analysing single sample data, set the resources to be used to filter variants with [gatk FilterVariantTranches](https://gatk.broadinstitute.org/hc/en-us/articles/360042479092-FilterVariantTranches`). For example:

```yaml
# Resources used to filter indels and SNP's
# Either call resources from databases defined above or provide a file directory
# Additional resources can be added
FILTERING:
  # Single sample analyses
  SINGLE: "--resource config['HAPMAP'])
          --resource config['MILLS']"
```

*Note: If not analysing single sample data, leave these fields blank*

If analysing cohort data, set the resources to be used to filter variants with [gatk VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042914791-VariantRecalibrator). For example:

```yaml
  # Cohort analyses
  COHORT:
    INDELS: "-resource:mills, known=false, training=true, truth=true, prior=12.0 config['MILLS'])
            -resource:dbsnp, known=true, training=false, truth=false, prior=2.0 config['dbSNP'])"
    SNPS: "-resource:hapmap, known=false,training=true,truth=true,prior=15.0 config['HAPMAP'])
          -resource:omni, known=false,training=true,truth=false,prior=12.0 config['OMNI'])
          -resource:1000G, known=false,training=true,truth=false,prior=10.0 config['SNP1000G'])
          -resource:dbsnp, known=true,training=false,truth=false,prior=2.0 config['dbSNP'])"
```

*Note: If not analysing cohort data, leave these fields blank*

Save your modified config file with a descriptive name

### 5. Create and activate a conda environment with python and snakemake and installed

```bash
conda create -n annot_pipeline_env python=3.7
conda activate annot_pipeline_env
conda install -c bioconda snakemake=5.14.0
```

### 6. Run the pipeline

Set the singularity bind location to a directory that contains the CADD database with the `--singularity-args` flag (eg. `'-B /home/lkemp/publicData/'`). Also specify your config file to be used with the `--configfile` flag and modify the number of cores to be used with the `-j` flag. First carry out a dry run. If there are no issues, start a full run without the `-n` flag.

Dry run:

```bash
snakemake -n -j 24 --use-conda --use-singularity --singularity-args '-B /dir/to/databases/' --configfile config.yaml
```

Full run:

```bash
snakemake -j 24 --use-conda --use-singularity --singularity-args '-B /dir/to/databases/' --configfile config.yaml
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html) for additional run parameters.

### 7. Evaluate the pipeline run

Generate an interactive html report

```bash
snakemake --report report.html --configfile config.yaml --report-stylesheet custom-stylesheet.css
```

### 8. Commit and push to your forked version of the repo

To maintain reproducibility, commit and push:

- All modified configuration file/s
- Output from the run such as reports and plots (optional)

### 9. Create a pull request with the [upstream repo](https://github.com/ESR-NZ/vcf_annotation_pipeline) to merge any useful changes (optional)

Contributions and feedback are more than welcome! :blush:

See [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) for help creating a pull request.

## Useful reading

Van der Auwera et al., (2013). *Current Protocols in Bioinformatics*. [From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/); 11(1110): 11.10.1–11.10.33.
