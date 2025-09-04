import os

from metasnek import fastq_finder, fasta_finder


# PARSE SAMPLES
config["koverage"]["samples"] = dict()
config["koverage"]["samples"]["reads"] = fastq_finder.parse_samples_to_dictionary(config["koverage"]["args"]["reads"])
config["koverage"]["samples"]["names"] = list(config["koverage"]["samples"]["reads"].keys())


# PARSE REF(S)
config["koverage"]["args"]["references"] = fasta_finder.parse_fastas(config["koverage"]["args"]["ref"])
if len(config["koverage"]["args"]["references"]) > 1:
    config["koverage"]["args"]["ref"] = os.path.join(config["koverage"]["args"]["temp"], "concatenated_refs.fasta")


# TARGETS
config["koverage"]["targets"] = dict()

if config["koverage"]["args"]["pafs"]:
    config["koverage"]["targets"]["pafs"] = expand(
        os.path.join(config["koverage"]["args"]["paf_dir"],"{sample}.paf.gz"),
        sample=config["koverage"]["samples"]["names"]
    )
else:
    config["koverage"]["targets"]["pafs"] = []

config["koverage"]["targets"]["coverage"] = [
    os.path.join(config["koverage"]["args"]["result"], "sample_coverage.tsv"),
    os.path.join(config["koverage"]["args"]["result"], "all_coverage.tsv"),
]

if config["koverage"]["args"]["report"]:
    config["koverage"]["targets"]["coverage"].append(os.path.join(config["koverage"]["args"]["result"], "report.html"))


config["koverage"]["targets"]["coverm"] = [
    os.path.join(config["koverage"]["args"]["result"], "sample_coverm_coverage.tsv")
]

config["koverage"]["targets"]["reports"] = [
    os.path.join(config["koverage"]["args"]["output"], "koverage.samples.tsv")
]


# KMER FILES
config["koverage"]["args"]["refkmers"] = os.path.join(
    config["koverage"]["args"]["temp"],
    os.path.basename(config["koverage"]["args"]["ref"]) + "." + str(config["koverage"]["args"]["kmer_size"]) + "mer.gz"
)
config["koverage"]["args"]["samplekmers"] = os.path.join(
    config["koverage"]["args"]["result"], "sample_kmer_coverage." + str(config["koverage"]["args"]["kmer_size"]) + "mer.tsv.gz"
)
config["koverage"]["args"]["allkmers"] = os.path.join(
    config["koverage"]["args"]["result"], "all_kmer_coverage." + str(config["koverage"]["args"]["kmer_size"]) + "mer.tsv.gz"
)

config["koverage"]["targets"]["kmercov"] = [
    config["koverage"]["args"]["samplekmers"],
    config["koverage"]["args"]["allkmers"]
]


# Add targets for pre-building the environments
config["koverage"]["targets"]["envs"] = []

for filename in os.listdir(os.path.join(workflow.basedir, "envs")):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        config["koverage"]["targets"]["envs"].append(os.path.join(config["koverage"]["args"]["temp"], filename + ".done"))
