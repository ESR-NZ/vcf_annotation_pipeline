vcf_annotation_pipeline annotates variant call format (VCF) data using `GATK4 <https://gatk.broadinstitute.org/hc/en-us>`_, `SnpSift <http://snpeff.sourceforge.net/SnpSift.html>`_, `VEP <https://asia.ensembl.org/info/docs/tools/vep/index.html>`_ and `genmod <http://moonso.github.io/genmod/>`_. 

Run parameters:
    * Reference genome build: {{ snakemake.config["BUILD"] }}
    * Input data type: {{ snakemake.config["DATA"] }} 
    * Sequencing type {{ snakemake.config["SEQUENCING"] }}