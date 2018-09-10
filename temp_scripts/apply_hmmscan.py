#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2018 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group

"""
This script:
1. Parse JSON-formatted parsed antiSMASH GBK output in specific filepath
  (or through direct API call)
2. Scan protein domains of each genes in a cluster, given that it is
  not done already (hmm_database and hmm_scanned properties) using
  specific <hmm_database.hmm>
3. Update the JSON-formatted file (or return the object) with hmmscan results

Example JSON output:
{
    "filename": "BGC0001.gbk",
    "hmm_database": "BiGFAM-merged-v0.0",
    "hmm_scanned": "05/05/2018 23:00:13",
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
                            "sequence": "AMAMAMAMAMAMAMA"
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
        print("apply_hmmscan.py <json_folder> <hmm_database>")
        sys.exit(1)

    json_folder = argv[1]
    hmm_database = argv[2]
        
    if os.path.isfile(json_folder):
        print("Error: {} is a file.".format(json_folder))
        sys.exit(1)

    json_files = []
    if os.path.isfile(json_folder):
        json_files = [json_folder]
    else:
        json_files = glob(os.path.join(json_folder,"**/*.json"), recursive=True)

    if not os.path.isfile(hmm_database):
        print("Warning: file {} did not exist.".format(hmm_database))
        sys.exit(1)

    print("Info: parsing HMM lengths (for alignment information)")
    hmm_len = get_hmm_lengths(hmm_database)

    for idx, json_file in enumerate(json_files):
        json_object = {}
        with open(json_file, "r") as json_text:
            try:                
                json_object = json.loads(json_text.read())
            except:
                print("Error: failed to load json file {}".format(json_file))
                continue
        result = apply_hmmscan(json_object, hmm_database, hmm_len, idx + 1)
        with open(json_file, "w") as json_text:
            json_text.write(json.dumps(result))

    end = time.time()
    print("Time elapsed: {}s".format(end - start))


def get_hmm_lengths(hmm_database_path):
    hmm_len = {}
    with open(hmm_database_path, "r") as hf:
        cur_name = ""
        for line in hf.readlines():
            if line.startswith("NAME"):
                cur_name = line.split(" ")[-1].rstrip()
            if line.startswith("LENG") and len(cur_name) > 0:
                hmm_len[cur_name] = int(line.split(" ")[-1].rstrip())
    return hmm_len


def apply_hmmscan(result, hmm_database_path, hmm_len, counter = 0):
    """
    given result JSON (proccessed or not processed),
    do hmmscan for every genes in the clusters using supplied hmm db
    """
    db_name = ".".join(hmm_database_path.split(".")[:-1]).split("/")[-1]

    #check previous run
    if ("hmm_database" in result) and (result["hmm_database"] == db_name):
        print("({}) Info: {} already hmmscanned on {}, skipping".format(counter, result["filename"], result["hmm_scanned"]))
        return result
        
    total_genes = 0
    for cluster in result["clusters"]:
        total_genes += len(cluster["genes"])
    if (total_genes < 1):
        print("({}) Error: {} have no genes, skipping".format(counter, result["filename"]))
        return result
    print("({}) Info: applying hmmscan on {} genes ({})".format(counter, total_genes, result["filename"]))

    # clear previous results
    for cluster in result["clusters"]:
        for gene in cluster["genes"]:
            gene["pfams"] = []

    concat_sequence = "" #concatenated sequences to speed up hmmsearch
    gene_starts = [] #start index of the genes in the concatenated sequences
    gene_records = [] #tuple of (cluster_idx, gene_idx) to trace back genes
    for cidx, cluster in enumerate(result["clusters"]):
        for gidx, gene in enumerate(cluster["genes"]):
            gene_starts.append(len(concat_sequence) + 1)
            gene_records.append((cidx, gidx))
            concat_sequence += gene["sequence"]
    hmm_results = hmmscan(concat_sequence, hmm_database_path, hmm_len)
    for hmm_result in hmm_results:
        start = hmm_result["seq_start"]
        end = hmm_result["seq_end"]
        for i, gene_start in enumerate(gene_starts):
            if (start > gene_start):
                if (i == len(gene_starts) - 1) or (end < gene_starts[i + 1]):
                    ci, gi = gene_records[i]
                    gene_record = result["clusters"][ci]["genes"][gi]
                    gene_record["pfams"].append({
                        "name": hmm_result["name"],
                        "bitscore": hmm_result["bitscore"],
                        "hit_start": (start - (gene_start - 1)),
                        "hit_end": (start - (gene_start - 1)) + (end - start),
                        "sequence": hmm_result["seq"]
                    });
                    break

    # apply timestamp
    result["hmm_database"] = db_name
    result["hmm_scanned"] = str(datetime.today())

    return result


def hmmscan(sequence, hmm_file, hmm_len):
    """
    given an AA sequence and a .hmm filepath, perform hmmscan
    and parse its result into JSON-formatted list
    example JSON:
    [
        {
            "model_start": 5,
            "model_end": 100,
            "seq_start": 10,
            "seq_end": 100,
            "bitscore": 230,
            "sequence": "AMAMAMAMAMMAMAMAA"
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
                for i in range(hsp.hit_end, hmm_len[hsp.hit_id] + 1):
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

    # filter overlapping hits
    temp_result = []
    result = sorted(result, key=lambda obj:obj["seq_start"])
    for hit in result:
        if len(temp_result) > 0:
            if (hit["seq_start"] - temp_result[-1]["seq_end"]) < 0: # is overlapping
                if hit["bitscore"] > temp_result[-1]["bitscore"]:
                    temp_result.pop()
                else:
                    continue
        temp_result.append(hit)
    result = temp_result
                
    return result


if __name__ == "__main__":
    main()