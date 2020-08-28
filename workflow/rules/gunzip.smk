rule gunzip:
    input:
        "../results/filtered/{sample}_filtered_multiallelicsites.vcf.gz"
    output:
        temp("../results/filtered/{sample}_filtered_multiallelicsites.vcf")
    conda:
        "../envs/gunzip.yaml"
    message:
        "Gunzipping {input}"
    shell:
        "gunzip < {input} > {output}"