# Run vcf_annotation_pipeline on a high performance cluster

## Table of contents

- [Run vcf_annotation_pipeline on a high performance cluster](#run-vcf_annotation_pipeline-on-a-high-performance-cluster)
  - [Table of contents](#table-of-contents)
  - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
  - [2. Take the pipeline to the data on the HPC](#2-take-the-pipeline-to-the-data-on-the-hpc)
  - [3. Setup files and directories](#3-setup-files-and-directories)
    - [Test data](#test-data)
  - [4. Get prerequisite software/hardware](#4-get-prerequisite-softwarehardware)
  - [5. Create a local copy of the GATK resource bundle (either b37 or hg38)](#5-create-a-local-copy-of-the-gatk-resource-bundle-either-b37-or-hg38)
    - [b37](#b37)
    - [hg38](#hg38)
  - [6. Create a local copy of other databases (either GRCh37 or GRCh38)](#6-create-a-local-copy-of-other-databases-either-grch37-or-grch38)
    - [GRCh37](#grch37)
    - [GRCh38](#grch38)
  - [7. Modify the configuration file](#7-modify-the-configuration-file)
    - [Overall workflow](#overall-workflow)
    - [Pipeline resources](#pipeline-resources)
    - [Variant filtering](#variant-filtering)
    - [VCF annotation](#vcf-annotation)
  - [8. Configure to run on a HPC](#8-configure-to-run-on-a-hpc)
  - [9. Modify the run scripts](#9-modify-the-run-scripts)
  - [10. Create and activate a conda environment with python and snakemake installed](#10-create-and-activate-a-conda-environment-with-python-and-snakemake-installed)
  - [11. Run the pipeline](#11-run-the-pipeline)
  - [12. Evaluate the pipeline run](#12-evaluate-the-pipeline-run)
  - [13. Commit and push to your forked version of the github repo](#13-commit-and-push-to-your-forked-version-of-the-github-repo)
  - [14. Repeat step 13 each time you re-run the analysis with different parameters](#14-repeat-step-13-each-time-you-re-run-the-analysis-with-different-parameters)
  - [15. Raise issues, create feature requests or create a pull request with the upstream repo to merge any useful changes to the pipeline (optional)](#15-raise-issues-create-feature-requests-or-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-to-the-pipeline-optional)

## 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help forking a repository

## 2. Take the pipeline to the data on the HPC

Clone the forked [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) repo into the same directory as your vcf files to be processed.

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help cloning a repository

## 3. Setup files and directories

*Note. for single sample analysis, you'll need the mapped bam files along with your vcf files. bam files are not needed for cohort analyses*

Required folder structure and file naming convention:

```bash

.
|___human_genomics_pipeline/results/called/
|                                      |___sample1_raw_snps_indels.vcf
|                                      |___sample1_raw_snps_indels.vcf.idx
|                                      |___sample2_raw_snps_indels.vcf
|                                      |___sample2_raw_snps_indels.vcf.idx
|                                      |___ ...
|
|___human_genomics_pipeline/results/mapped/
|                                      |___sample1_recalibrated.bam
|                                      |___sample1_recalibrated.bai
|                                      |___sample2_recalibrated.bam
|                                      |___sample2_recalibrated.bai
|                                      |___ ...
|
|___vcf_annotation_pipeline/

```

If you're analysing cohort's of samples, you will need a directory with a [pedigree file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) for each cohort/family using the following folder structure and file naming convention:

```bash

.
|___human_genomics_pipeline/results/called/
|                                      |___sample1_raw_snps_indels.g.vcf
|                                      |___sample1_raw_snps_indels.g.vcf.idx
|                                      |___sample2_raw_snps_indels.g.vcf
|                                      |___sample2_raw_snps_indels.g.vcf.idx
|                                      |___sample3_raw_snps_indels.g.vcf
|                                      |___sample3_raw_snps_indels.g.vcf.idx
|                                      |___sample4_raw_snps_indels.g.vcf
|                                      |___sample4_raw_snps_indels.g.vcf.idx
|                                      |___sample5_raw_snps_indels.g.vcf
|                                      |___sample5_raw_snps_indels.g.vcf.idx
|                                      |___sample6_raw_snps_indels.g.vcf
|                                      |___sample6_raw_snps_indels.g.vcf.idx
|                                      |___ ...
|
|___pedigrees/
|     |___proband1_pedigree.ped
|     |___proband2_pedigree.ped
|     |___ ...
|
|___vcf_annotation_pipeline/

```

Requirements:

- Currently, the filenames of the pedigree files need to be labelled with the name of the proband/individual affected with the disease phenotype in the cohort (we will be working towards removing this requirement)
- Singletons and cohorts need to be run in separate pipeline runs
- It is assumed that there is one proband/individual affected with the disease phenotype of interest in a given cohort (one individual with a value of 2 in the 6th column of a given pedigree file)

### Test data

The provided [test dataset](./test) can be used. Setup the test dataset before running the pipeline on this data - choose to setup to run either a single sample analysis or a cohort analysis with the `-a` flag. For example:

```bash
cd ./human_genomics_pipeline
bash ./test/setup_test.sh -a cohort
```

## 4. Get prerequisite software/hardware

For GPU accelerated runs, you'll need [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) and [NVIDIA CLARA PARABRICKS and dependencies](https://www.nvidia.com/en-us/docs/parabricks/local-installation/). Talk to your system administrator to see if the HPC has this hardware and software available.

Other software required to get setup and run the pipeline:

- [Git](https://git-scm.com/) (tested with version 2.7.4)
- [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) (tested with version 4.8.2)
- [Mamba](https://github.com/TheSnakePit/mamba) (tested with version 0.4.4) (note. [mamba can be installed via conda with a single command](https://mamba.readthedocs.io/en/latest/installation.html#existing-conda-install))
- [gsutil](https://pypi.org/project/gsutil/) (tested with version 4.52)
- [gunzip](https://linux.die.net/man/1/gunzip) (tested with version 1.6)

Most of this software is commonly pre-installed on HPC's, likely available as modules that can be loaded. Talk to your system administrator if you need help with this.

## 5. Create a local copy of the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) (either b37 or hg38)

### b37

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37?prefix=)

```bash
gsutil cp -r gs://gatk-legacy-bundles/b37 /where/to/download/
```

### hg38

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

```bash
gsutil cp -r gs://genomics-public-data/resources/broad/hg38 /where/to/download/
```

## 6. Create a local copy of other databases (either GRCh37 or GRCh38)

### GRCh37

Download the [Ensembl-VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) database using a [conda version of Ensembl-VEP](https://anaconda.org/bioconda/ensembl-vep)

```bash
conda create -n download_data_env python=3.7
conda activate download_data_env
conda install -c bioconda ensembl-vep=99.2
vep_install -a cf -s homo_sapiens -y GRCh37 -c /output/file/path/GRCh37 --CONVERT
conda deactivate
```

Download the [CADD database](https://cadd.gs.washington.edu/download) and it's associated index file.

```bash
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
```

Create a custom [dbNSFP database](https://sites.google.com/site/jpopgen/dbNSFP) build by following [this documentation](https://github.com/GenomicsAotearoa/dbNSFP_build)

### GRCh38

Download [Ensembl-VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) database using a [conda install of Ensembl-VEP](https://anaconda.org/bioconda/ensembl-vep)

```bash
mamba create -n download_data_env python=3.7
conda activate download_data_env
mamba install -c bioconda ensembl-vep=99.2
vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/file/path/GRCh38 --CONVERT
conda deactivate
```

Download the [CADD database](https://cadd.gs.washington.edu/download) and it's associated index file.

```bash
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

Create a custom [dbNSFP database](https://sites.google.com/site/jpopgen/dbNSFP) build by following [this documentation](https://github.com/GenomicsAotearoa/dbNSFP_build)

## 7. Modify the configuration file

Edit 'config.yaml' found within the config directory.

### Overall workflow

Specify the build of reference genome (either 'GRCh37' (for b37) or 'GRCh38' (for hg38)) to use

```yaml
BUILD: "GRCh37"
```

Specify whether the data is to be analysed on it's own ('Single') or as a part of a cohort of samples ('Cohort'). For example:

```yaml
DATA: "Single"
```

Specify whether the pipeline should be GPU accelerated where possible (either 'Yes' or 'No', this requires [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) and [NVIDIA CLARA PARABRICKS](https://www.nvidia.com/en-us/docs/parabricks/local-installation/))

```yaml
GPU_ACCELERATED: "No"
```

Specify whether the data should be prepared to be ingested into [scout](http://www.clinicalgenomics.se/scout/) (either 'Yes' or 'No'). This will output an additional vcf file that has undergone additional filtering steps. For example:

```yaml
PREPARE_FOR_SCOUT: "Yes"
```

Set the working directories to the reference human genome file (b37 or hg38). For example:

```yaml
REFGENOME: "/scratch/publicData/b37/human_g1k_v37_decoy.fasta"
```

Set the the working directory to your dbSNP database file (b37 or hg38). For example:

```yaml
dbSNP: "/scratch/publicData/b37/dbsnp_138.b37.vcf"
```

Set the the working directory to a temporary file directory. Make sure this is a location with a fair amount of memory space for large intermediate analysis files. For example:

```yaml
TEMPDIR: "/scratch/tmp/"
```

If analysing WES data, pass a design file (.bed) indicating the genomic regions that were sequenced (see [here](https://leahkemp.github.io/documentation/human_genomic_pipelines/design_files.html) for more information on accessing design files). Also set the level of padding by passing the amount of padding in base pairs. For example:

*If NOT analysing WES data, leave these fields blank*

```yaml
WES:
  # File path to the exome capture regions over which to operate
  INTERVALS: "/scratch/publicData/sure_select_human_all_exon_V7/S31285117_Padded.bed"
  # Padding (in bp) to add to each region
  PADDING: "100"
```

### Pipeline resources

These settings allow you to configure the resources per rule/sample

Set the number of threads to use per sample/rule for multithreaded rules (`rule gatk_CNNScoreVariants`, `rule vep` and `genmod models`). Multithreading will significantly speed up these rules, however the improvements in speed will diminish beyond 8 threads. If desired, a different number of threads can be set for these multithreaded rules by utilising the `--set-threads` flag in the runscript (see [step 8](#8-modify-the-run-scripts)).

```yaml
THREADS: 8
```

Set the maximum memory usage per rule/sample (eg. '40g' for 40 gigabytes, this should suffice for exomes)

```yaml
MAXMEMORY: "40g"
```

Set the maximum number of GPU's to be used per rule/sample for gpu-accelerated runs (eg `1` for 1 GPU)

```yaml
GPU: 1
```

It is a good idea to consider the number of samples that you are processing. For example, if you set `THREADS: "8"` and set the maximum number of cores to be used by the pipeline in the run script to `-j 32` (see step 6), a maximum of 3 samples will be able to run at one time for these rules (if they are deployed at the same time), but each sample will complete faster. In contrast, if you set `THREADS: "1"` and `-j 32`, a maximum of 32 samples could be run at one time, but each sample will take longer to complete. This also needs to be considered when setting `MAXMEMORY` + `--resources mem_mb` and `GPU` + `--resources gpu`.

### Variant filtering

If analysing single sample data, pass the resources to be used to filter variants with [gatk FilterVariantTranches](https://gatk.broadinstitute.org/hc/en-us/articles/360042479092-FilterVariantTranches`) to the `--resource` flag. For example:

*If NOT analysing single sample data, leave these fields blank*

```yaml
FILTERING:
  # ...for analysis of single samples
  SINGLE: "--resource /scratch/publicData/b37/hapmap_3.3.b37.vcf
          --resource /scratch/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
```

If analysing cohort data, pass the resources to be used to filter variants with [gatk VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042914791-VariantRecalibrator) to the `--resource` flag. For example:

*If NOT analysing cohort data, leave these fields blank*

```yaml
  # ...for analysis of cohorts
  COHORT:
    INDELS: "--resource:mills,known=false,training=true,truth=true,prior=12.0 /scratch/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 /scratch/publicData/b37/1000G_phase1.indels.b37.vcf
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /scratch/publicData/b37/dbsnp_138.b37.vcf"
    SNPS: "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /scratch/publicData/b37/hapmap_3.3.b37.vcf
          --resource:omni,known=false,training=true,truth=false,prior=12.0 /scratch/publicData/b37/1000G_omni2.5.b37.vcf
          --resource:1000G,known=false,training=true,truth=false,prior=10.0 /scratch/publicData/b37/1000G_phase1.indels.b37.vcf
          --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /scratch/publicData/b37/dbsnp_138.b37.vcf"
```

If running as GPU accelerated, use the following syntax when setting the resources to be used to filter variants with [gatk VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042914791-VariantRecalibrator)

```yaml
  COHORT:
    INDELS: "--resource mills,known=false,training=true,truth=true,prior=12.0:/scratch/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
            --resource 1000G,known=false,training=true,truth=false,prior=10.0:/scratch/publicData/b37/1000G_phase1.indels.b37.vcf
            --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/scratch/publicData/b37/dbsnp_138.b37.vcf"
    SNPS: "--resource hapmap,known=false,training=true,truth=true,prior=15.0:/scratch/publicData/b37/hapmap_3.3.b37.vcf
          --resource omni,known=false,training=true,truth=false,prior=12.0:/scratch/publicData/b37/1000G_omni2.5.b37.vcf
          --resource 1000G,known=false,training=true,truth=false,prior=10.0:/scratch/publicData/b37/1000G_phase1.indels.b37.vcf
          --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/scratch/publicData/b37/dbsnp_138.b37.vcf"
```

Set the tranche filtering level for snps and indels (by [gatk FilterVariantTranches](https://gatk.broadinstitute.org/hc/en-us/articles/360041417412-FilterVariantTranches) for single samples and [gatk VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360041851391-VariantRecalibrator) for cohorts). For example:

```yaml
  TRANCHE: 
    SNPS: "99.95"
    INDELS: "99.4"
```

### VCF annotation

Set the the working directories to the other vcf annotation databases (GRCh37 or GRCh38). For example:

```yaml
VEP: "/scratch/publicData/vep/GRCh37/"
dbNSFP: "/scratch/publicData/dbNSFP/GRCh37/dbNSFPv4.0a.hg19.custombuild.gz"
CADD: "/scratch/publicData/CADD/GRCh37/whole_genome_SNVs.tsv.gz"
```

## 8. Configure to run on a HPC

*This will deploy the non-GPU accelerated rules to slurm and deploy the GPU accelerated rules locally (pbrun_cnnscorevariants). Therefore, if running the pipeline gpu accelerated, the pipeline should be deployed from the machine with the GPU's.*

In theory, this cluster configuration should be adaptable to other job scheduler systems. In fact, snakemake is designed to allow pipelines/workflows such as this to be portable to different job schedulers by having any [any cluster specific configured outside of the workflow/pipeline](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated). but here I will provide a minimal example for deploying this pipeline to [slurm](https://slurm.schedmd.com/).

Configure `account:` and `partition:` in the default section of 'cluster.json' in order to set the parameters for slurm sbatch (see documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated) and [here](https://slurm.schedmd.com/)). For example:

```json
{
    "__default__" :
    {
        "account" : "lkemp",
        "partition" : "prod",
        "output" : "logs/slurm-%j_{rule}_{wildcards.sample}.out"
    }
}
```

There are a plethora of additional slurm parameters that can be configured (and can be configured per rule). If you set additional slurm parameters, remember to pass them to the `--cluster` flag in the runscripts. See [here](https://snakemake-on-nesi.sschmeier.com/snake.html#slurm-and-nesi-specific-setup) for a good working example of deploying a snakemake workflow to [NeSi](https://www.nesi.org.nz/).

## 9. Modify the run scripts

Set the singularity bind location to a directory that contains your pipeline working directory with the `--singularity-args '-B'` flag. Set the number maximum number of cores to be used with the `--cores` flag and the maximum amount of memory to be used (in megabytes) with the `resources mem_mb=` flag. If running GPU accelerated, also set the maximum number of GPU's to be used with the `--resources gpu=` flag. For example:

Dry run (dryrun_hpc.sh):

```bash
snakemake \
--dryrun \
--cores 32 \
--resources mem_mb=150000 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--singularity-args '-B /scratch/' \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
-o {cluster.output}"
```

Full run (run_hpc.sh):

```bash
snakemake \
--cores 32 \
--resources mem_mb=150000 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--singularity-args '-B /scratch/' \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
-o {cluster.output}"
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html#all-options) for additional run parameters.

## 10. Create and activate a conda environment with python and snakemake installed

```bash
cd ./workflow/
mamba env create -f pipeline_run_env.yml
conda activate pipeline_run_env
```

## 11. Run the pipeline

First carry out a dry run

```bash
bash dryrun_hpc.sh
```

If there are no issues, start a full run

```bash
bash run_hpc.sh
```

## 12. Evaluate the pipeline run

Generate an interactive html report

```bash
bash report.sh
```

## 13. Commit and push to your forked version of the github repo

To maintain reproducibility, commit and push:

- All configuration files
- All run scripts
- The final report

## 14. Repeat step 13 each time you re-run the analysis with different parameters

## 15. Raise issues, create feature requests or create a pull request with the [upstream repo](https://github.com/ESR-NZ/vcf_annotation_pipeline) to merge any useful changes to the pipeline (optional)

See [the README](https://github.com/ESR-NZ/vcf_annotation_pipeline/blob/dev/README.md#contribute-back) for info on how to contribute back to the pipeline!
