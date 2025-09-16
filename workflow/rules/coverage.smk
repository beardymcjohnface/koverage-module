import os

rule koverage_idx_ref:
    """Prepare the reference fasta file for mapping with minimap"""
    input:
        config["koverage"]["args"]["ref"]
    output:
        config["koverage"]["args"]["ref"] + '.idx'
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        **config["resources"]["med"]
    conda:
        os.path.join("..", "envs", "minimap.yaml")
    container:
        config["koverage"]["container"]["minimap"]
    envmodules:
        *config["koverage"]["envmodules"]["minimap"]
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "idx_ref.txt")
    log:
        os.path.join(config["koverage"]["args"]["output_paths"]["log"], "idx_ref.err")
    shell:
        ("awk 'BEGIN {{count=-1}} /^>/ {{ $0 = \">\" ++count }} 1' {input} "
            "| minimap2 -t {threads} -d {output} - 2> {log}")


rule koverage_faidx_ref:
    """Index the reference fasta file with samtools faidx"""
    input:
        config["koverage"]["args"]["ref"]
    output:
        config["koverage"]["args"]["ref"] + '.fai'
    conda:
        os.path.join("..", "envs", "minimap.yaml")
    container:
        config["koverage"]["container"]["minimap"]
    envmodules:
        *config["koverage"]["envmodules"]["minimap"]
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "faidx_ref.txt")
    log:
        os.path.join(config["koverage"]["args"]["output_paths"]["log"], "faidx_ref.err")
    shell:
        "samtools faidx {input} 2> {log}"


rule koverage_raw_coverage:
    """Map and collect the raw read counts for each sample
    
    lib: single line, total number of reads
    counts: "contig\tcontig_len\tcount\tmean\tmedian\thitrate\tvariance
    """
    input:
        ref = config["koverage"]["args"]["ref"] + ".idx",
        r1=lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R1"],
        fai = config["koverage"]["args"]["ref"] + '.fai'
    output:
        counts = temp(os.path.join(config["koverage"]["args"]["output_paths"]["temp"], "{sample}.counts.pkl")),
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        **config["resources"]["med"]
    params:
        r2 = lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] if config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] else "",
        pafs = config["koverage"]["args"]["pafs"],
        paf_dir = config["koverage"]["args"]["output_paths"]["paf"],
        bin_width = config["koverage"]["args"]["bin_width"],
        minimap = config["koverage"]["params"]["minimap"]
    conda:
        os.path.join("..", "envs", "minimap.yaml")
    container:
        config["koverage"]["container"]["minimap"]
    envmodules:
        *config["koverage"]["envmodules"]["minimap"]
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "raw_coverage.{sample}.txt")
    log:
        err = os.path.join(config["koverage"]["args"]["output_paths"]["log"], "raw_coverage.{sample}.err"),
    script:
        os.path.join("..", "scripts", "minimapWrapper.py")


rule koverage_sample_coverage:
    """convert raw counts to coverage values
    
    output: sample\tcontig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance
    """
    input:
        counts = os.path.join(config["koverage"]["args"]["output_paths"]["temp"],"{sample}.counts.pkl"),
    output:
        temp(os.path.join(config["koverage"]["args"]["output_paths"]["temp"],"{sample}.cov.tsv"))
    params:
        binwidth = config["koverage"]["args"]["bin_width"]
    threads: 1
    conda:
        os.path.join("..", "envs", "numpy.yaml")
    container:
        config["koverage"]["container"]["numpy"]
    envmodules:
        *config["koverage"]["envmodules"]["numpy"]
    log:
        err =os.path.join(config["koverage"]["args"]["output_paths"]["log"], "sample_coverage.{sample}.err"),
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "sample_coverage.{sample}.txt")
    script:
        os.path.join("..", "scripts", "sampleCoverage.py")


rule koverage_all_sample_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(os.path.join(config["koverage"]["args"]["output_paths"]["temp"],"{sample}.cov.tsv"), sample=config["koverage"]["samples"]["names"])
    output:
        os.path.join(config["koverage"]["args"]["output_paths"]["results"], "sample_coverage.tsv")
    threads: 1
    log:
        err = os.path.join(config["koverage"]["args"]["output_paths"]["log"], "all_sample_coverage.err"),
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "all_sample_coverage.txt")
    shell:
        ("printf 'Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance\n' > {output} 2> {log}; "
        "cat {input} >> {output} 2> {log} ")


rule koverage_combine_coverage:
    """Combine all sample coverages"""
    input:
        coverage = os.path.join(config["koverage"]["args"]["output_paths"]["results"],"sample_coverage.tsv"),
        fai = config["koverage"]["args"]["ref"] + '.fai'
    output:
        all_cov = os.path.join(config["koverage"]["args"]["output_paths"]["results"], "all_coverage.tsv"),
    threads: 1
    log:
        err = os.path.join(config["koverage"]["args"]["output_paths"]["log"], "combine_coverage.err"),
    benchmark:
        os.path.join(config["koverage"]["args"]["output_paths"]["bench"], "combine_coverage.txt")
    script:
        os.path.join("..", "scripts", "combineCoverage.py")
