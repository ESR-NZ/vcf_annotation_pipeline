if config['DATA'] == "Single" or config['DATA'] == 'single':
    infile = "../results/filtered/{sample}_filtered_scoutfiltered.vcf.gz"
elif config['DATA'] == "Cohort" or config['DATA'] == 'cohort':
    infile = "../results/filtered/{sample}_filtered_multiallelicsites.vcf.gz"

if config['DATA'] == "Single" or config['DATA'] == 'single':
    outfile = "../results/filtered/{sample}_filtered_scoutfiltered.vcf"
elif config['DATA'] == "Cohort" or config['DATA'] == 'cohort':
    outfile = "../results/filtered/{sample}_filtered_multiallelicsites.vcf.gz"

rule gunzip:
    input:
        infile
    output:
        temp(outfile)
    conda:
        "../envs/gunzip.yaml"
    message:
        "Gunzipping {input}"
    shell:
        "gunzip < {input} > {output}"