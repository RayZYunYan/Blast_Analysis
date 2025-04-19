from Bio.Blast import NCBIWWW
from Bio import SeqIO
import time
import os

fasta_files = [f for f in os.listdir("source_seq") if f.endswith(".fasta")]
lengths = [100, 500, 1000]

for fasta in fasta_files:
    record = SeqIO.read(f"source_seq/{fasta}", "fasta")
    seq_id = os.path.splitext(fasta)[0]

    for L in lengths:
        sequence = str(record.seq[:L])
        label = f"{seq_id}_{L}bp"

        for algo, db in [("blastn", "nt"), ("blastx", "nr")]:
            print(f"\nRunning {algo.upper()} for {label}...")
            start = time.time()
            result_handle = NCBIWWW.qblast(algo, db, sequence)
            end = time.time()

            os.makedirs("output_seq", exist_ok=True)
            filename = f"output_seq/{algo}_{label}.xml"
            with open(filename, "w") as out:
                out.write(result_handle.read())

            print(f"Saved: {filename} | Time: {end - start:.2f}s")
