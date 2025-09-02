import os


rule koverage_sample_tsv:
    output:
        tsv = os.path.join(config["koverage"]["args"]["output"], "koverage.samples.tsv")
    params:
        sample_dict = config["koverage"]["samples"]["reads"]
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params.sample_dict, output.tsv)


rule koverage_report:
    """Generate html report for koverage"""
    input:
        smpl = os.path.join(config["koverage"]["args"]["result"], "sample_coverage.tsv"),
        all = os.path.join(config["koverage"]["args"]["result"], "all_coverage.tsv")
    output:
        html = os.path.join(config["koverage"]["args"]["result"], "report.html")
    params:
        sample_cov_desc = config["koverage"]["report"]["map"]["sample_cov_desc"],
        all_cov_desc = config["koverage"]["report"]["map"]["all_cov_desc"],
        sample_names = config["koverage"]["samples"]["names"],
        ref_fasta = config["koverage"]["args"]["ref"],
        max_ctg = config["koverage"]["args"]["report_max_ctg"]
    threads: 1
    conda:
        os.path.join("..", "envs", "report.yaml")
    container:
        config["koverage"]["container"]["report"]
    envmodules:
        *config["koverage"]["envmodules"]["report"]
    log:
        err = os.path.join(config["koverage"]["args"]["log"], "coverage_report.err"),
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"], "coverage_report.txt")
    script:
        os.path.join("..", "scripts", "koverageReport.py")


rule koverage_build_env:
    output:
        os.path.join(config["koverage"]["args"]["temp"], "{env}.done")
    conda:
        lambda wildcards: os.path.join("..", "envs", wildcards.env)
    shell:
        "touch {output}"
