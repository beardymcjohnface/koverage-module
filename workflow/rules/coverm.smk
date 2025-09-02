import os


rule koverage_coverm_map_pe:
    input:
        ref = config["koverage"]["args"]["ref"],
        r1=lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R1"],
    output:
        bam = os.path.join(config["koverage"]["args"]["temp"], "{sample}.bam"),
        bai = os.path.join(config["koverage"]["args"]["temp"], "{sample}.bam.bai")
    params:
        r2 = lambda wildcards: config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] if config["koverage"]["samples"]["reads"][wildcards.sample]["R2"] else "",
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
        os.path.join(config["koverage"]["args"]["bench"], "coverm_map_pe.{sample}.txt")
    log:
        os.path.join(config["koverage"]["args"]["log"], "coverm_map_pe.{sample}.err")
    shell:
        "{{ "
        "minimap2 "
            "-t {threads} "
            "-ax sr "
            "--secondary=no "
            "{input.ref} "
            "{input.r1} "
            "{params.r2} "
        "| samtools sort "
            "-T {wildcards.sample} "
            "-@ {threads} - "
        "| samtools view "
            "-bh -F 4 "
            "> {output.bam} ;"
        "samtools index "
            "{output.bam}; "
        "}} 2> {log}"


rule koverage_coverm_bam2counts:
    input:
        os.path.join(config["koverage"]["args"]["temp"], "{sample}.bam")
    output:
        os.path.join(config["koverage"]["args"]["temp"], "{sample}.cov")
    params:
        params = config["koverage"]["params"]["coverm"]
    conda:
        os.path.join("..", "envs", "coverm.yaml")
    container:
        config["koverage"]["container"]["coverm"]
    envmodules:
        *config["koverage"]["envmodules"]["coverm"]
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"],"coverm_bam2counts.{sample}.txt")
    log:
        os.path.join(config["koverage"]["args"]["log"], "coverm_bam2counts.{sample}.err")
    shell:
        "coverm contig "
            "-b {input} " 
            "{params.params} "
            "> {output} "
            "2> {log} "


rule koverage_coverm_combine:
    input:
        expand(os.path.join(config["koverage"]["args"]["temp"], "{sample}.cov"), sample=config["koverage"]["samples"]["names"])
    output:
        os.path.join(config["koverage"]["args"]["result"], "sample_coverm_coverage.tsv")
    params:
        samples = config["koverage"]["samples"]["names"],
        dir = config["koverage"]["args"]["temp"]
    run:
        with open(output[0], "w") as outfh:
            with open(input[0], "r") as infh:
                header = infh.readline().split("\t")
                header = [' '.join(element.split()[1:]) for element in header]
                header[0] = "Contig"
                header = '\t'.join(header)
            outfh.write("Sample\t" + header + "\n")
            for sample in params.samples:
                with open(os.path.join(params.dir, sample + ".cov"), "r") as infh:
                    infh.readline()
                    for line in infh:
                        outfh.write(sample + "\t" + line)


rule koverage_reneo_coverage:
    input:
        expand(os.path.join(config["koverage"]["args"]["temp"],"{sample}.cov"),sample=config["koverage"]["samples"]["names"])
    output:
        os.path.join(config["koverage"]["args"]["result"], "reneo.coverage.tsv")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        **config["resources"]["ram"]
    benchmark:
        os.path.join(config["koverage"]["args"]["bench"],"reneo_coverage.txt")
    log:
        os.path.join(config["koverage"]["args"]["log"], "reneo_coverage.err")
    shell:
        ("for i in {input}; do "
            "tail -n+2 $i; "
        "done | "
            "awk -F '\t' '{{ sum[$1] += $2 }} END {{ for (key in sum) print key, sum[key] }}' "
            "> {output} 2> {log}; ")
