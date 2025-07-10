#!/usr/bin/env python3

import os
import subprocess
import csv
from Bio.Blast import NCBIXML

# Paths
assemblies_dir = "assemblies"
query_file = "pezAT_refs.fasta"

# list of assemblies
assemblies = [f for f in os.listdir(assemblies_dir) if f.endswith(".fa")]

# prepare summary list
summary = []

for genome in assemblies:
    genome_path = os.path.join(assemblies_dir, genome)
    db_name = genome_path.replace(".fa", "_db")

    # make blast db
    subprocess.run(["makeblastdb", "-in", genome_path, "-dbtype", "nucl", "-out", db_name],
                   check=True)

    # run blastn
    output_xml = genome_path.replace(".fa", "_pezAT.xml")
    blastn_cmd = [
        "blastn",
        "-query", query_file,
        "-db", db_name,
        "-outfmt", "5",
        "-out", output_xml
    ]
    subprocess.run(blastn_cmd, check=True)

    # parse blast xml
    with open(output_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        # initialize counters
        pezA_count = 0
        pezT_count = 0

        print(f"\nResults for {genome}:")

        for blast_record in blast_records:
            # robustly identify gene name in query
            query = blast_record.query.lower()
            if "peza" in query:
                gene = "pezA"
            elif "pezt" in query:
                gene = "pezT"
            else:
                gene = "unknown"

            if blast_record.alignments:
                for align in blast_record.alignments:
                    for hsp in align.hsps:
                        identity_pct = (hsp.identities / hsp.align_length) * 100
                        print(f"  Query: {gene}")
                        print(f"    Hit: {align.hit_def}")
                        print(f"    % identity: {identity_pct:.1f}")
                        print(f"    alignment length: {hsp.align_length}")
                        print(f"    subject coords: {hsp.sbjct_start}-{hsp.sbjct_end}")

                        if gene == "pezA":
                            pezA_count += 1
                        elif gene == "pezT":
                            pezT_count += 1
            else:
                print(f"  Query: {gene}")
                print("    No hits found.")

        # add to summary
        summary.append({
            "genome": genome,
            "pezA_hits": pezA_count,
            "pezT_hits": pezT_count
        })

# write summary table
with open("pezAT_hits_summary.csv", "w", newline="") as csvfile:
    fieldnames = ["genome", "pezA_hits", "pezT_hits"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in summary:
        writer.writerow(row)

print("\nSummary table saved to pezAT_hits_summary.csv.")
print("BLAST pezAT check complete.")
