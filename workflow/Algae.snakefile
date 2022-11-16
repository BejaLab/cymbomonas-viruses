
rule algae_getorf:
    input:
        "algae/{genome}.fna"
    output:
        "analysis/algae/{genome}.faa"
    params:
        minsize = 300,
        maxsize = 100000
    conda:
        "envs/emboss.yaml"
    shell:
        "getorf -minsize {params.minsize} -maxsize {params.maxsize} -filter {input} > {output}"
