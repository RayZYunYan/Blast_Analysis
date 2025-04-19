from Bio.Blast import NCBIXML
import pandas as pd
import glob
import re
import os

records = []

for xml_file in glob.glob("output_seq/blast[nx]_*.xml"):
    match = re.match(r"output_seq/(blast[nx])_(\w+)_([0-9]+)bp\.xml", xml_file.replace("\\", "/"))
    if not match:
        continue

    algorithm, filename, length = match.groups()

    with open(xml_file) as handle:
        try:
            blast_record = NCBIXML.read(handle)
        except Exception:
            continue

        if not blast_record.alignments:
            continue

        top_hit = blast_record.alignments[0]
        hsp = top_hit.hsps[0]

        evalue = hsp.expect
        identity = 100 * hsp.identities / hsp.align_length
        score = hsp.score

        records.append({
            "algorithm": algorithm,
            "file": filename,
            "length": int(length),
            "evalue": evalue,
            "identity_percent": round(identity, 2),
            "score": score
        })

df = pd.DataFrame(records)
df = df.sort_values(by=["file", "length", "algorithm"])
df.to_csv("blast_summary.csv", index=False)
print(df)
