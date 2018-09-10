#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2018 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group

"""
This script:
1. Parse JSON-formatted, hmmscanned antiSMASH GBK output in specific filepath
  (or through direct API call)
2. Given a list of core pfam domains, accompanied by their respective subdomain
hmm models (in a folder), scan every core domain hits sequence in a BGC
3. Update the JSON-formatted file (or return the object) with hmmscan results of
the subdomains

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
- HMMER ver...
- BioPython ver...
- helperlibs ver...
"""


import sys
import os
import json
from datetime import datetime
from glob import glob
from sys import argv
from Bio import SeqIO
from Bio import SearchIO
from helperlibs.wrappers.io import TemporaryDirectory
import subprocess


def main():
    import time
    start = time.time()

    if len(argv) < 2:
        print("apply_subdomainscan.py <json_folder> <coredomains_list> <subdomain_database_folder>")
        sys.exit(1)

    json_folder = argv[1]
    coredomains_list_txt = argv[2]
    subdomain_database_folder = argv[3]
        
    if os.path.isfile(json_folder):
        print("Error: {} is a file.".format(json_folder))
        sys.exit(1)

    json_files = []
    if os.path.isfile(json_folder):
        json_files = [json_folder]
    else:
        json_files = glob(os.path.join(json_folder,"**/*.json"), recursive=True)

    if not os.path.isfile(coredomains_list_txt):
        print("Warning: file {} did not exist.".format(coredomains_list_txt))
        sys.exit(1)
        
    if os.path.isfile(json_folder):
        print("Error: {} is a file.".format(json_folder))
        sys.exit(1)

    coredomains_list = []
    with open(coredomains_list_txt, "r") as cf:
        for line in cf.readlines():
            coredomains_list.append(line.rstrip())

    for idx, json_file in enumerate(json_files):
        json_object = {}
        with open(json_file, "r") as json_text:
            try:                
                json_object = json.loads(json_text.read())
            except:
                print("Error: failed to load json file {}".format(json_file))
                continue
        result = apply_subdomainscan(json_object, coredomains_list, subdomain_database_folder, idx + 1)
        with open(json_file, "w") as json_text:
            json_text.write(json.dumps(result))

    end = time.time()
    print("Time elapsed: {}s".format(end - start))


def apply_subdomainscan(result, coredomains_list, subdomain_database_folder, counter = 0):
    """
    given result JSON (proccessed or not processed),
    do hmmscan for every genes in the clusters using supplied hmm db
    """
    db_name = subdomain_database_folder.split("/")[-1]

    #check previous run
    if ("subdomain_database" in result) and (result["subdomain_database"] == db_name):
        print("({}) Info: {} already subdomainscanned on {}, skipping".format(counter, result["filename"], result["subdomain_scanned"]))
        return result

    #check if hmmscanned
    if "hmm_database" not in result:
        print("({}) Info: {} haven't got hmmscanned, skipping".format(counter, result["filename"]))
        return result

    print("({}) Info: applying subdomainscan on {}".format(counter, result["filename"]))

    # clear previous results and apply subdomainscan for respective hmm hits
    for cluster in result["clusters"]:
        for gene in cluster["genes"]:
            if "pfams" in gene:
                for pfam in gene["pfams"]:
                    if "subdomains" in pfam:
                        del(pfam["subdomains"])
                    if pfam["name"] in coredomains_list:
                        subdomain_hmm = os.path.join(subdomain_database_folder, "{}-subdomains.hmm".format(pfam["name"]))
                        pfam["subdomains"] = []
                        for res in subdomainscan(pfam["sequence"], subdomain_hmm, len(pfam["sequence"])):
                            pfam["subdomains"].append({
                                "name": res["name"],
                                "start": res["seq_start"],
                                "end": res["seq_end"],
                                "bitscore": res["bitscore"],
                                "sequence": res["seq"]
                            })

    # apply timestamp
    result["subdomain_database"] = db_name
    result["subdomain_scanned"] = str(datetime.today())

    return result


def subdomainscan(sequence, hmm_file, hmm_len):
    """
    given an AA sequence and a .hmm filepath, perform hmmscan
    and parse its results into JSON-formatted list
    example JSON:
    [
        {
            "name": "Pfam0001",
            "model_start": 5,
            "model_end": 100,
            "seq_start": 10,
            "seq_end": 100,
            "bitscore": 230,
            "seq": "AMAMAMAMAMMAMAMAA"
        }
    ]
    """
    result = []

    with TemporaryDirectory() as temp_dir:
        fasta_path = os.path.join(temp_dir, "temp.fa")
        fasta_file = open(fasta_path, "w")
        fasta_file.write(">temp\n{}\n".format(sequence))
        fasta_file.close()
        result_path = os.path.join(temp_dir, "temp.output")
        command = "hmmscan -o {} -T 20 --incT 20 {} {}".format(result_path, hmm_file, fasta_path)
        subprocess.check_output(command, shell=True)
        for runresult in SearchIO.parse(result_path, 'hmmer3-text'):
            for hsp in runresult.hsps:
                padding_left = ""
                padding_right = ""
                for i in range(0, hsp.hit_start):
                    padding_left += "-"
                for i in range(hsp.hit_end, hmm_len + 1):
                    padding_right += "-"
                if hsp.hit_strand != hsp.query_strand:
                    padding_left, padding_right = padding_right, padding_left
                seq = ""
                hmmseq = str(hsp.hit.seq)
                queryseq = str(hsp.query.seq)
                for i in range(0, len(hmmseq)):
                    if hmmseq[i] != '.':
                        seq += queryseq[i]
                seq = "".join([padding_left, seq, padding_right])
                result.append({
                    "name": hsp.hit_id,
                    "model_start": hsp.hit_start,
                    "model_end": hsp.hit_end,
                    "seq_start": hsp.query_start,
                    "seq_end": hsp.query_end,
                    "bitscore": hsp.bitscore,
                    "seq": seq
                })

    # filter hits (of the same subdomain)
    temp_result = []
    result = sorted(result, key=lambda obj:obj["bitscore"], reverse=True)
    added_hits = []
    for hit in result:
        if hit["name"] in added_hits:
            continue
        else:
            added_hits.append(hit["name"])
            temp_result.append(hit)
    result = temp_result
                
    return result


if __name__ == "__main__":
    main()