#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2018 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group
"""
This script:
1. Parse antiSMASH GBK outputs in specific folder
2. Convert them to JSON-formatted files ready to be processed by other scripts

Example JSON output:
{
    "filename": "BGC0001.gbk",
    "clusters": [
        {
            "name": "cluster001",
            "class": "",
            "genes": [
                {
                    "name": "gene001",
                    "start": 1,
                    "end": 61,
                    "sequence": "AMAMAMAMAMAMAMAMMAMAMAMAMAMAMAMAMAMAMMAMAMAMAMAMAMAMAMAMAMMA"
                },
                ...
            ]
        },
        ...
    ]
}

Requirements:
- BioPython ver...
"""


import sys
import os
import json
from glob import glob
from sys import argv
from Bio import SeqIO


def main():
    """
    arguments: build_json.py <gbks_folder> <output_folder>
    """
    gbks_folder = argv[1]
    output_folder = argv[2]

    if os.path.isfile(output_folder):
        print("Error: {} is a file.".format(output_folder))
        sys.exit(1)

    for result in load_gbks(gbks_folder):
        base_filename = ".".join(result["filename"].split(".")[:-1]).split("/")[-1]
        json_filepath = os.path.join(output_folder, "{}.json".format(base_filename))
        with open(json_filepath, "w") as jsonfile:
            jsonfile.write(json.dumps(result))


def load_gbks(path):
    """
    load all gbks in a folder (or load one gbk file),
    return list of JSON-formatted BGCs
    example json output:
    { see script description }
    """
    results = []
    filenames = []

    if os.path.isfile(path):
        filenames = [path]
    else:
        filenames = glob(os.path.join(path,"**/*.gbk"), recursive=True)

    for filename in filenames:
        records = []
        try:
            records = list(SeqIO.parse(filename, "genbank"))
            if len(records) > 1:
                print("Info: file {} contains more than 1 record".format(filename))
        except:
            print("Error: file {} is not a correct gbk file.".format(filename))
            continue
            # sys.exit(1)
        for record in records:
            clusters = fetch_clusters(record)
            results.append({
                "filename": filename.split("/")[-1],
                "clusters": [{
                    "genes": clusters
                }]
            })

    return results


def fetch_clusters(record):
    """
    given a biopython SeqRecord, fetch all clusters information
    and return list of JSON-formatted BGCs
    example json output:
    { see script description }
    - todo: extract already present domain information from antiSMASH (for efficiency)
    """
    clusters = []
    cds_list = []

    for feature in record.features:
        if feature.type in ["cluster"]:
            clusters.append({
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "class": [] if "product" not in feature.qualifiers else feature.qualifiers["product"][0].split("-")
            })
        elif feature.type == "CDS":
            cds = {
                "info": {},
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "sequence": feature.qualifiers["translation"][0]
            }
            included_info = ["locus_tag", "gene", "protein_id"]
            for qualifier in included_info:
                if qualifier in feature.qualifiers:
                    cds["info"][qualifier] = feature.qualifiers[qualifier][0]
            cds_list.append(cds)
    
    return cds_list

    if len(clusters) > 1:
        print("Info: record {} contains more than 1 cluster".format(record.id))

    cds_list = sorted(cds_list, key = lambda cds: cds["start"], reverse = True)

    for cluster in clusters:
        cluster_cds = []
        for cds in cds_list:
            if (cds["start"] < cluster["start"]) or (cds["end"] > cluster["end"]):
                break
            cluster_cds.append(cds)
        cluster["genes"] = cluster_cds

    return clusters

if __name__ == "__main__":
    main()