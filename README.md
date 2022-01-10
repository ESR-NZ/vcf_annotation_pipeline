# vcf_annotation_pipeline

A Snakemake workflow to filter raw variants (snp and indels) and annotate vcf (variant call format) files (single samples or cohorts) using [GATK4](https://gatk.broadinstitute.org/hc/en-us), [SnpSift](http://snpeff.sourceforge.net/SnpSift.html), [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html), [genmod](https://github.com/moonso/genmod) and [dbSNP](https://www.ncbi.nlm.nih.gov/SNP/). The vcf file can also optionally be prepared for ingestion into [scout](http://www.clinicalgenomics.se/scout/) which involves some filtering steps and/or run on [NVIDIA GPU's](https://www.nvidia.com/en-gb/graphics-cards/) where [nvidia clara parabricks software is available](https://www.nvidia.com/en-us/docs/parabricks/quickstart-guide/software-overview/) for speedups in analysis times. This pipeline is designed to follow [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and the data ingested into [scout](https://github.com/Clinical-Genomics/scout) for clinical interpretation. However, this pipeline also stands on it's own, as a vcf annotation pipeline. This pipeline has been developed with human genetic data in mind, however we designed it to be species agnostic. Genetic data from other species can be analysed by setting a species-specific reference genome and variant databases in the configuration file (but not all situations have been tested).

- [vcf_annotation_pipeline](#vcf_annotation_pipeline)
  - [Pipeline summary - single samples](#pipeline-summary---single-samples)
  - [Pipeline summary - single samples - GPU accelerated](#pipeline-summary---single-samples---gpu-accelerated)
  - [Pipeline summary - cohort samples](#pipeline-summary---cohort-samples)
  - [Pipeline summary - cohort samples - GPU accelerated](#pipeline-summary---cohort-samples---gpu-accelerated)
  - [Main output files](#main-output-files)
  - [Prerequisites](#prerequisites)

## Pipeline summary - single samples

<img src="./images/rulegraph_single.png" class="center">

## Pipeline summary - single samples - GPU accelerated

<img src="./images/rulegraph_single_gpu.png" class="center">

## Pipeline summary - cohort samples

<img src="./images/rulegraph_cohort.png" class="center">

## Pipeline summary - cohort samples - GPU accelerated

<img src="./images/rulegraph_cohort_gpu.png" class="center">

## Main output files

Single samples:

- `results/filtered/sample1_filtered.vcf`
- `results/annotated/sample1_filtered_annotated.vcf`
- `results/readyforscout/sample1_filtered_annotated_readyforscout.vcf.gz`

Cohort samples:

- `results/filtered/sample1_filtered.vcf`
- `results/annotated/sample1_filtered_annotated.vcf`
- `results/readyforscout/sample1_filtered_annotated_readyforscout.vcf`

## Prerequisites

- **Prerequisite hardware:** [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) (for GPU accelerated runs)
- **Prerequisite software:** [NVIDIA CLARA parabricks and dependencies](https://www.nvidia.com/en-us/docs/parabricks/local-installation/) (for GPU accelerated runs), [Git](https://git-scm.com/) (tested with version 2.7.4), [Mamba](https://github.com/TheSnakePit/mamba) (tested with version 0.4.4) with [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) (tested with version 4.8.2), [gsutil](https://pypi.org/project/gsutil/) (tested with version 4.52), [gunzip](https://linux.die.net/man/1/gunzip) (tested with version 1.6), [R](https://www.r-project.org/) (tested with version 3.2.2)
