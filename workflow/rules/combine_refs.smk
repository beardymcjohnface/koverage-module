import os
from metasnek import fasta_finder


rule koverage_combine_fastas:
    """Combine the multiple ref fastas"""
    output:
        os.path.join(config["koverage"]["args"]["temp"],"concatenated_refs.fasta")
    params:
        ref_files = config["koverage"]["args"]["references"]
    localrule:
        True
    run:
        fasta_finder.combine_fastas(params.ref_files, output[0])
