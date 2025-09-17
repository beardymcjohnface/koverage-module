import os

from metasnek import fastq_finder

# PARSE SAMPLES
config["koverage"]["samples"] = dict()
config["koverage"]["samples"]["reads"] = fastq_finder.parse_samples_to_dictionary(config["koverage"]["args"]["reads"])
config["koverage"]["samples"]["names"] = list(config["koverage"]["samples"]["reads"].keys())
