
rule extract_segments:
    input:
        "output/viral_segments.gbk"

rule extract_element:
    input:
        lambda w: expand("analysis/segments/{host}/segments.gbk", host = segments[w.segment]['genome'])
    output:
        "analysis/PLV_segments/{segment}.gbk"
    params:
        flank = segment_flanks,
        segment = lambda w: segments[w.segment]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/extract_element.py"

rule mcl_add_to_gbk:
    input:
        data = "output/mcl_genes.tsv",
        gbk = expand("analysis/PLV_segments/{segment}.gbk", segment = segment_names)
    output:
        "output/viral_segments.gbk"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/mcl_add_to_gbk.py"

