from Bio.Blast import NCBIWWW
from Bio import SeqIO
import time
import os
from datetime import datetime

# Dynamically generated log files
log_time = datetime.now().strftime("%Y-%b-%d-%H-%M-%S")
log_path = f"Blast_Analysis/logs/{log_time}.log"
os.makedirs("Blast_Analysis/logs", exist_ok=True)

def write_log_entry(label, algo, status, runtime):
    with open(log_path, "a") as log:
        log.write(f"Running content logs: algorithms: {algo.upper()}，Objective: {label}\n")
        log.write(f"Operational status: {status}\n")
        log.write(f"Runtime: {runtime if runtime != '-' else '-'}s\n\n")

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

            try:
                start = time.time()
                result_handle = NCBIWWW.qblast(algo, db, sequence)
                end = time.time()

                os.makedirs("output_seq", exist_ok=True)
                filename = f"output_seq/{algo}_{label}.xml"
                with open(filename, "w") as out:
                    out.write(result_handle.read())

                elapsed = round(end - start, 2)
                print(f"Saved: {filename} | Time: {elapsed}s")
                write_log_entry(label, algo, "✅ Done", elapsed)

            except Exception as e:
                print(f"❌ Error running {algo.upper()} for {label}: {e}")
                write_log_entry(label, algo, "❌ Not Done", "-")
