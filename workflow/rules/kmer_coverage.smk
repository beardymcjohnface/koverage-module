import os


rule koverage_jellyfish_db:
    """Calculate a jellyfish database of the reads"""
    input:
        r1=lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R1"],
    output:
        os.path.join(config["koverage"]["args"]["temp"], "{sample}." + str(config["koverage"]["args"]["kmer_size"]) + "mer"),
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        **config["resources"]["med"]
    params:
        r2 = lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] if config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] else "",
        kmer = config["koverage"]["args"]["kmer_size"],
        jf = config["koverage"]["params"]["jellyfish"],
        cat = lambda wildcards: "gunzip -c" if config["koverage"]["samples"]["reads"][wildcards.sample]["R1"].endswith(".gz") else "cat"
    conda:
        os.path.join("..", "envs", "jellyfish.yaml")
    envmodules:
        *config["koverage"]["envmodules"]["jellyfish"]
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "jellyfish_db.{sample}.txt")
    log:
        os.path.join(config["koverage"]["args"]["log"], "jellyfish.{sample}.err")
    shell:
        ("jellyfish count "
            "-m {params.kmer} "
            "-t {threads} "
            "-o {output} "
            "{params.jf} "
            "<({params.cat} {input.r1} {params.r2}) "
            "2> {log}")


rule koverage_ref_kmer_prep:
    """Sample kmers for a reference fasta"""
    input:
        config["koverage"]["args"]["ref"]
    output:
        config["koverage"]["args"]["refkmers"]
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        **config["resources"]["med"]
    params:
        ksize = config["koverage"]["args"]["kmer_size"],
        kspace = config["koverage"]["args"]["kmer_sample"],
        kmin = config["koverage"]["args"]["kmer_min"],
        kmax = config["koverage"]["args"]["kmer_max"],
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "ref_kmer_prep.txt")
    log:
        err = os.path.join(config["koverage"]["args"]["log"], "ref_kmer_prep.err"),
    script:
        os.path.join("..", "scripts", "refSampleKmer.py")


rule koverage_kmer_screen:
    """Screen jellyfish database for ref kmers"""
    input:
        ref = config["koverage"]["args"]["refkmers"],
        db = os.path.join(config["koverage"]["args"]["temp"], "{sample}." + str(config["koverage"]["args"]["kmer_size"]) + "mer")
    output:
        temp(os.path.join(config["koverage"]["args"]["temp"], "{sample}." + str(config["koverage"]["args"]["kmer_size"]) + "mer.kcov.gz"))
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        **config["resources"]["med"]
    conda:
        os.path.join("..", "envs","jellyfish.yaml")
    envmodules:
        *config["koverage"]["envmodules"]["jellyfish"]
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "kmer_screen.{sample}.txt")
    log:
        err = os.path.join(config["koverage"]["args"]["log"], "kmer_screen.{sample}.err"),
    script:
        os.path.join("..", "scripts", "kmerScreen.py")


rule koverage_all_sample_kmer_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(
            os.path.join(
                config["koverage"]["args"]["temp"],
                "{sample}." + str(config["koverage"]["args"]["kmer_size"]) + "mer.kcov.gz"
            ),
            sample=config["koverage"]["samples"]["names"]
        )
    output:
        config["koverage"]["args"]["samplekmers"]
    threads: 1
    # conda:
    #     os.path.join("..", "envs", "zstd.yaml")
    # envmodules:
    #     *config["koverage"]["envmodules"]["zstd"]
    log:
        os.path.join(config["koverage"]["args"]["log"], "all_sample_kmer_coverage.err")
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "all_sample_kmer_coverage.txt")
    shell:
        ("{{ "
            "printf 'Sample\tContig\tSum\tMean\tMedian\tHitrate\tVariance\n' 2> {log}; "
            "zcat {input} 2> {log}; "
        "}} | gzip -1 - > {output} ")


rule koverage_combine_kmer_coverage:
    """Combine all sample kmer coverages"""
    input:
        config["koverage"]["args"]["samplekmers"]
    output:
        all_cov = config["koverage"]["args"]["allkmers"]
    threads: 1
    log:
        err = os.path.join(config["koverage"]["args"]["log"], "combine_kmer_coverage.err"),
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "combine_kmer_coverage.txt")
    script:
        os.path.join("..", "scripts", "combineKmerCoverage.py")
