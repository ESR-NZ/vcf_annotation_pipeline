rule gunzip:
    input:
        "../results/readyforscout/{sample}_filtered_annotated_multiallelicsites.vcf.gz"
    output:
        temp("../results/readyforscout/{sample}_filtered_annotated_multiallelicsites.vcf")
    conda:
        "../envs/gunzip.yaml"
    message:
        "Gunzipping {input}"
    shell:
        "gunzip < {input} > {output}"