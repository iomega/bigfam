#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2018 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group
"""
This script:
1. Parse JSON-formatted, hmmscanned, subdomainscanned antiSMASH GBK output
in specific filepath (or through direct API call)
2. extract a feature matrix consisted of domains presence + subdomains presence
and update the JSON files

Example JSON output:
{
    "filename": "BGC0001.gbk",
    "hmm_database": "BiGFAM-merged-v0.0",
    "hmm_scanned": "05/05/2018 23:00:13",
    "subdomain_database": "BiGFAM-merged-v0.0-subdomains",
    "subdomain_scanned": "05/06/2018 17:20:53",
    "clusters": [
        {
            "name": "cluster001",
            "features": [1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0.5, 0.23, 1, 0. 0.15, 0.78]
            "genes": [
                {
                    "name": "gene001",
                    "start": 1,
                    "end": 61,
                    "sequence": "AMAMAMAMAMAMAMAMMAMAMAMAMAMAMAMAMAMAMMAMAMAMAMAMAMAMAMAMAMMA"
                    "pfams": [
                        {
                            "name": "PF001",
                            "start": 1,
                            "end": 16,
                            "bitscore": 123,
                            "sequence": "AMAMAMAMAMAMAMA",
                            "subdomains": [
                                {
                                    "name": "PF001-c1",
                                    "start": 1,
                                    "end": 10,
                                    "bitscore": 200,
                                    "sequence": "AMAMAMAMAM"
                                }
                            ]
                        }
                    ]
                },
                ...
            ]
        },
        ...
    ]
}

Requirements:
"""


import sys
import os
import json
from datetime import datetime
from glob import glob
from sys import argv



def main():
    """
    arguments: calculate_features.py <json_folder> <hmm_database_path> <subdomain_database_folder>
    """
    json_folder = argv[1]
    hmm_database_path = argv[2]
    subdomain_database_folder = argv[3]

    # load jsons
    json_files = []
    if os.path.isfile(json_folder):
        json_files = [json_folder]
    else:
        json_files = glob(os.path.join(json_folder,"**/*.json"), recursive=True)

    # fetch domain names
    domain_list = []
    with open(hmm_database_path, "r") as hmm_file:
        for line in hmm_file.readlines():
            if line.startswith("NAME"):
                domain_list.append(line.split(" ")[-1].rstrip())
    domain_list = sorted(list(set(domain_list)))

    # fetch subdomain names
    subdomain_list = []
    subdomain_files = []
    if os.path.isfile(subdomain_database_folder):
        subdomain_files = [subdomain_database_folder]
    else:
        subdomain_files = glob(os.path.join(subdomain_database_folder,"**/*.hmm"), recursive=True)
    for subdomain_file in subdomain_files:
        with open(subdomain_file, "r") as hmm_file:
            for line in hmm_file.readlines():
                if line.startswith("NAME"):
                    subdomain_list.append(line.split(" ")[-1].rstrip())
    subdomain_list = sorted(list(set(subdomain_list)))

    # extract features
    for json_file in json_files:
        json_object = {}
        with open(json_file, "r") as json_text:
            try:                
                json_object = json.loads(json_text.read())
            except:
                print("Error: failed to load json file {}".format(json_file))
                continue
        result = extract_features(json_object, domain_list, subdomain_list)
        with open(json_file, "w") as json_text:
            json_text.write(json.dumps(result))


def extract_features(result, domain_list, subdomain_list):
    result["features"] = []
    print("Info: extracting features ({})".format(result["filename"]))
    result["features"].extend(count_domains(result, domain_list))
    result["features"].extend(count_subdomains(result, subdomain_list))
    return result


def count_domains(result, domain_list):
    scores = [0] * len(domain_list)
    for cluster in result["clusters"]:
        for gene in cluster["genes"]:
            if "pfams" in gene:
                for hit in gene["pfams"]:
                    scores[domain_list.index(hit["name"])] += 1
    maximum = max(1, max(scores))
    scores = [float(score / maximum) for score in scores]
    return scores


def count_subdomains(result, subdomain_list):
    scores_bgc = [0.00] * len(subdomain_list)
    for cluster in result["clusters"]:
        for gene in cluster["genes"]:
            scores_per_domain = {}
            if "pfams" in gene:
                for hit in gene["pfams"]:
                    if "subdomains" in hit:
                        if hit["name"] not in scores_per_domain:
                            scores_per_domain[hit["name"]] = [0] * len(subdomain_list)
                        for subhit in hit["subdomains"]:
                            scores_per_domain[hit["name"]][subdomain_list.index(subhit["name"])] = max(scores_per_domain[hit["name"]][subdomain_list.index(subhit["name"])], subhit["bitscore"])
            scores = [0.00] * len(subdomain_list)
            for domain in scores_per_domain:
                maximum = max(1, max(scores_per_domain[domain]))
                scores_per_domain[domain] = [float(score / maximum) for score in scores_per_domain[domain]]
                for i in range(0, len(subdomain_list)):
                    scores[i] = max(scores[i], scores_per_domain[domain][i])
            for i in range(0, len(subdomain_list)):
                scores_bgc[i] = max(scores_bgc[i], scores[i])
    return scores_bgc


if __name__ == "__main__":
    main()