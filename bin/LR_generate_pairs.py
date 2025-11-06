#!/usr/bin/env python3

import os
import re
import gzip
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import edlib
import re
from collections import namedtuple, Counter
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description="Detect and change reverse linker seq")
parser.add_argument("--fastq",help="fastq input")
parser.add_argument("--new_fname",help="New filename")
parser.add_argument("--out_dir",help="Output directory")
parser.add_argument("--l1dist",type=int,help="Edit distance for first linker")
parser.add_argument("--l2dist",type=int,help="Edit distance for second linker")
parser.add_argument("--chemistry",help="Kit chemistry version (e.g. v1, v2, v3)")
parser.add_argument("--multiple_fq",action="store_true",help="Multiple fastq files")
parser.add_argument("--max_deletions",default=1,type=int,help="Max deletions for edit distance")
parser.add_argument("--threads",default=1,type=int,help="Number of threads available.")

args = parser.parse_args()

chem = args.chemistry
out_dir = args.out_dir
filename = args.fastq
l1dist = args.l1dist
l2dist = args.l2dist
new_fname = args.new_fname
multiple = args.multiple_fq
max_d = args.max_deletions

if chem == "v1" or chem == "v2":
    L1F = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
    L2F = "ATCCACGTGCTTGAGACTGTGG"
    L1RC = "AGTCGTACGCCGATGCGAAACATCGGCCAC"
    L2RC = "CCACAGTCTCAAGCACGTGGAT"
elif chem == "v3":
    L1F = "ATGAGGGGTCAG"
    L2F = "TCCAACCACCTC"
    L1RC = "CTGACCCCTCAT"
    L2RC = "GAGGTGGTTGGA"
else:
    print("Chemistry not set")  # Assuming you want to print this message

r1fwd_list=[]
r2fwd_list=[]
r1rev_list=[]
r2rev_list=[]

if multiple:
    set_num = re.findall("_[0-9]*\.", filename)[0]
else:
    set_num = "."

#cnts = {"n_reads":0, "rc_match":0,"fwd_match":0, "no_match":0,
#    "no_bc3":0, "rev_bad_pos":0, "fwd_bad_pos":0}

def clean_header(record):
    record.description = record.description.split(" ")[0]

Match = namedtuple('Match', ['start', 'end', 'dist'])

def _count_deletions(cigar):
    return sum(int(n) if n else 1 for n, op in re.findall(r'(\d*)([MID])', cigar) if op == 'D')

def fnm_edlib(query, sequence, max_l_dist, max_deletions):
    matches = []

    result = edlib.align(query, sequence, mode='HW', task='locations', k=max_l_dist)
    if result['editDistance'] == -1:
        return []

    for start, end in result.get('locations', []):
        subseq = sequence[start:end+1]
        detailed = edlib.align(query, subseq, mode='NW', task='path')

        #Extract cigar string:
        cigar_string = detailed.get('cigar', '')

        #Only proceed if cigar is not None. None means that the alignment failed
        if cigar_string is not None:
            deletions = _count_deletions(detailed.get('cigar', ''))
            if deletions <= max_deletions:
                matches.append(Match(start=start, end=end+1, dist=detailed['editDistance']))

    return matches


#IMPORTANT NOTE: this won't work in the case there are two matches with the minimal distance as one will be choosen randomly
def extract_closest_match(matches):
    # Step 1: Find the minimum distance
    min_dist = min(m.dist for m in matches)

    # Step 2: Get all matches with that minimum distance
    best_matches = [m for m in matches if m.dist == min_dist]

    return best_matches[0]

def process_record(record, L1RC, L2RC, L1F, L2F, l1dist, l2dist, max_d):
    result = {
        "r1": None,
        "r2": None,
        "orientation": None,
        "cnts": {"rc_match": 0, "fwd_match": 0, "rev_bad_pos": 0, "fwd_bad_pos": 0, "no_bc3": 0, "no_match": 0},
        "unmatched": None
    }

    #Extract the sequence:
    sequence = record.seq

    #Set control variable (boolean):
    accepted = False

    #Search nearest matches for linkers sequences:
    #Reverse complement configuration:
    l1rc_match = fnm_edlib(L1RC, sequence, l1dist, max_d)

    if l1rc_match:
        #Get the best match:
        best_l1rc = extract_closest_match(l1rc_match)

        #Window to subset:
        start_subset = best_l1rc.start - 30
        end_subset = best_l1rc.start

        #Only proceed to L2RC if L1RC found:
        l2rc_match = fnm_edlib(L2RC, sequence[start_subset:end_subset], l2dist, max_d)

    else:
        l2rc_match = []

    if len(l1rc_match) and len(l2rc_match) > 0:

        #Get the best match:
        best_l2rc = extract_closest_match(l2rc_match)

        #Check if linkers are in the proper configuration:
        if best_l1rc.start - (best_l2rc.end + start_subset) == 8:
            clean_header(record)
            end = best_l1rc.end + 18
            r1rev = record[1:(end - 57)]
            r2rev = record[(end - 58):end]

            if len(r2rev) < 58:
                result["cnts"]["rev_bad_pos"] += 1
            elif end <= len(record):
                r2rev.seq = r2rev.seq.reverse_complement()
                result.update({"r1": r1rev, "r2": r2rev, "orientation": "rc"})
                result["cnts"]["rc_match"] += 1
            else:
                result["cnts"]["no_bc3"] += 1

            #Update control variable:
            accepted = True

    #Forward configuration:
    if not accepted:
        l1f_match = fnm_edlib(L1F, sequence, l1dist, max_d)

        if l1f_match:
            #Get the best match:
            best_l1f = extract_closest_match(l1f_match)

            #Window to subset:
            start_subset = best_l1f.end
            end_subset = best_l1f.end + 30

            #Only proceed to L2F if L1F found:
            l2f_match = fnm_edlib(L2F, sequence[start_subset:end_subset], l2dist, max_d)

        else:
            l2f_match = []

        if len(l1f_match) and len(l2f_match) > 0:

            #Get the best match:
            best_l2f = extract_closest_match(l2f_match)

            if (best_l2f.start + start_subset) - best_l1f.end == 8:
                clean_header(record)
                start = best_l1f.start - 18
                r1fwd = record[(start + 58):len(record)]
                r2fwd = record[start:(start + 58)]

                if len(r2fwd) < 58:
                    result["cnts"]["fwd_bad_pos"] += 1
                elif start > 0:
                    r1fwd.seq = r1fwd.seq.reverse_complement()
                    result.update({"r1": r1fwd, "r2": r2fwd, "orientation": "fwd"})
                    result["cnts"]["fwd_match"] += 1
                else:
                    result["cnts"]["no_bc3"] += 1

                #Update control variable:
                accepted = True

    if not accepted:
        result["cnts"]["no_match"] += 1
        #SeqIO.write(record, out_dir + "unmatched.fastq", "fastq")
        result["unmatched"] = record

    return result

if not os.path.exists(filename):
    print(f"File: {e} doesn't exist, skipping")
else:
    try:
        with open(filename, "rt") as handle:

            records = SeqIO.parse(handle, "fastq")

            for record in records:

                results = Parallel(n_jobs=args.threads, prefer="processes")(
                delayed(process_record)(record, L1RC, L2RC, L1F, L2F, l1dist, l2dist, max_d)
                for record in records
                )

            #Post-process results:
            r1rev_list, r2rev_list = [], []
            r1fwd_list, r2fwd_list = [], []
            cnts = Counter()
            unmatched = []

            for res in results:
                cnts.update(res["cnts"])
                if res["orientation"] == "rc":
                    r1rev_list.append(res["r1"])
                    r2rev_list.append(res["r2"])
                elif res["orientation"] == "fwd":
                    r1fwd_list.append(res["r1"])
                    r2fwd_list.append(res["r2"])
                if res["unmatched"]:
                    unmatched.append(res["unmatched"])


            #Write PE-like reads:
            r1_complete = r1fwd_list + r1rev_list
            r2_complete = r2fwd_list + r2rev_list
            out_r1 = out_dir + f"{new_fname}_R1" + set_num + "fastq.gz"
            out_r2 = out_dir + f"{new_fname}_R2" + set_num + "fastq.gz"
            out_unmatched = out_dir + f"{new_fname}_unmatched.fastq.gz"

            #Create output files:
            with gzip.open(out_r1, "wt") as h1, gzip.open(out_r2, "wt") as h2, gzip.open(out_unmatched, "wt") as h3:

                SeqIO.write(r1_complete, h1, "fastq")
                SeqIO.write(r2_complete, h2, "fastq")
                SeqIO.write(unmatched, h3, "fastq")

            print(filename, out_r1, out_r2)
            print(cnts)
            print("Fraction of reads retained",round((cnts["rc_match"]+cnts["fwd_match"])/sum(cnts.values()),3), "\n")

    except Exception as exptn:
        print(f"An error occurred while processing the file: {exptn}")
