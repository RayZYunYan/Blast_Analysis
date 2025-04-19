# 🔬 BLAST Identity Comparison Tool

This project compares the identity performance of BLASTn and BLASTx algorithms across DNA sequences of varying lengths and complexity.

---

## 📁 Project Structure

```
BLAST_ANALYSIS/
├── source_seq/             # Input .fasta sequences
│   ├── 1sequence.fasta
│   ├── 2sequence.fasta
│   └── ...
├── output_seq/             # All BLAST XML results (auto-generated)
│   ├── blastn_1sequence_100bp.xml
│   └── ...
├── venv/                   # Python virtual environment
├── blast_summary.csv       # Final summary table (auto-generated)
├── batch_run_blast_multi.py     # Main BLAST script
├── summarize_results.py         # Parses XML to generate CSV
├── parse_blast_results.py       # Optional: view top hits in terminal
├── plot_results.py              # Generates 3 comparison graphs
├── plot_length_vs_identity.png          # Output graph 1
├── plot_algorithm_identity_avg.png      # Output graph 2
├── plot_grouped_bar_length_algo.png     # Output graph 3
```

---

## 💻 Setup Instructions

### 1. Create and activate virtual environment (Windows)

```bash
python -m venv venv
.\env\Scripts\activate
```

### 2. Install required libraries

```bash
pip install biopython pandas seaborn matplotlib
```

---

## 🚀 Run the Workflow

### Step 1: Perform BLASTn and BLASTx queries

```bash
python batch_run_blast_multi.py
```

This reads `source_seq/*.fasta`, runs BLASTn/BLASTx at 100/500/1000bp, and saves results to `output_seq/`.

---

### Step 2: Summarize results into CSV

```bash
python summarize_results.py
```

This parses all `.xml` results in `output_seq/` and outputs `blast_summary.csv`.

---

### Step 3: Generate visualizations

```bash
python plot_results.py
```

Creates 3 graph files:
- `plot_length_vs_identity.png`
- `plot_algorithm_identity_avg.png`
- `plot_grouped_bar_length_algo.png`

---

## 📊 Sample Results Overview

- **BLASTn** achieves consistently high identity across all sequence lengths.
- **BLASTx** shows more variability, especially at 500bp.
- Longer sequences do not always improve match quality for BLASTx.

---

## 📌 Conclusion

This tool helps explore when BLASTn or BLASTx performs better depending on sequence length and complexity. It can be easily extended to more sequences or to include BLASTp.

---

## ✍️ Author

CS123A Term Project  
Spring 2025

Ray Zhang: rui.zhang@sjsu.edu  
Jisheng Jiang: jisheng.jiang@sjsu.edu
