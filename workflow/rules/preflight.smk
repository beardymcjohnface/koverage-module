import os

from metasnek import fasta_finder


"""
Expand output dir and file paths to include base output dir
"""
for output_file_path in config["koverage"]["args"]["output_paths"]:
    config["koverage"]["args"]["output_paths"][output_file_path] =  os.path.join(
        config["koverage"]["args"]["output"], config["koverage"]["args"]["output_paths"][output_file_path]
    )


"""
PARSE REF(S)
"""
config["koverage"]["args"]["references"] = fasta_finder.parse_fastas(config["koverage"]["args"]["ref"])
if len(config["koverage"]["args"]["references"]) > 1:
    config["koverage"]["args"]["ref"] = os.path.join(config["koverage"]["args"]["temp"], "concatenated_refs.fasta")


"""
TARGETS
"""
config["koverage"]["targets"] = dict()

if config["koverage"]["args"]["pafs"]:
    config["koverage"]["targets"]["pafs"] = expand(
        os.path.join(config["koverage"]["args"]["output_paths"]["paf"],"{sample}.paf.gz"),
        sample=config["koverage"]["samples"]["names"]
    )
else:
    config["koverage"]["targets"]["pafs"] = []

config["koverage"]["targets"]["coverage"] = [
    os.path.join(config["koverage"]["args"]["output_paths"]["results"], "sample_coverage.tsv"),
    os.path.join(config["koverage"]["args"]["output_paths"]["results"], "all_coverage.tsv"),
]

if config["koverage"]["args"]["report"]:
    config["koverage"]["targets"]["coverage"].append(os.path.join(config["koverage"]["args"]["output_paths"]["results"], "report.html"))


config["koverage"]["targets"]["coverm"] = [
    os.path.join(config["koverage"]["args"]["output_paths"]["results"], "sample_coverm_coverage.tsv")
]

config["koverage"]["targets"]["reports"] = [
    os.path.join(config["koverage"]["args"]["output"], "koverage.samples.tsv")
]


"""
KMER FILES
"""
config["koverage"]["args"]["refkmers"] = os.path.join(
    config["koverage"]["args"]["output_paths"]["temp"],
    os.path.basename(config["koverage"]["args"]["ref"]) + "." + str(config["koverage"]["args"]["kmer_size"]) + "mer.gz"
)
config["koverage"]["args"]["samplekmers"] = os.path.join(
    config["koverage"]["args"]["output_paths"]["results"], "sample_kmer_coverage." + str(config["koverage"]["args"]["kmer_size"]) + "mer.tsv.gz"
)
config["koverage"]["args"]["allkmers"] = os.path.join(
    config["koverage"]["args"]["output_paths"]["results"], "all_kmer_coverage." + str(config["koverage"]["args"]["kmer_size"]) + "mer.tsv.gz"
)

config["koverage"]["targets"]["kmercov"] = [
    config["koverage"]["args"]["samplekmers"],
    config["koverage"]["args"]["allkmers"]
]


"""
Add targets for pre-building the environments - TODO: delete if migrating to snk
"""
config["koverage"]["targets"]["envs"] = []

for filename in os.listdir(os.path.join(workflow.basedir, "envs")):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        config["koverage"]["targets"]["envs"].append(os.path.join(config["koverage"]["args"]["output_paths"]["temp"], filename + ".done"))
