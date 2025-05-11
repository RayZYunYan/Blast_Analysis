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
├── pics/                   # Visualization output directory
│   ├── plot_algorithm_identity_avg.png
│   ├── plot_grouped_bar_length_algo.png
│   ├── plot_identity_compare.png
│   └── plot_length_vs_identity.png
├── venv/                   # Python virtual environment
├── blast_summary.csv       # Final summary table (auto-generated)
├── batch_run_blast_multi.py     # Main BLAST script
├── summarize_results.py         # Parses XML to generate CSV
├── parse_blast_results.py       # Optional: view top hits in terminal
├── plot_results.py              # Generates comparison graphs
├── BLAST_Colab.ipynb            # Google Colab notebook version
```

---

## 💻 Setup Instructions

### 1. Create and activate virtual environment (Windows)

```bash
python -m venv venv
.\venv\Scripts\activate
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

Creates graph files saved in the `/pics` directory:
- `plot_length_vs_identity.png`
- `plot_algorithm_identity_avg.png`
- `plot_grouped_bar_length_algo.png`

For the charts, we generate them locally via log files as we have a total of 2 devices, which reduces the production time. The images for this part are stored in the `/pics` directory.

---

## 📊 Sample Results Overview

- **BLASTn** achieves consistently high identity across all sequence lengths.
- **BLASTx** shows more variability, especially at 500bp.
- Longer sequences do not always improve match quality for BLASTx.

---

## 💻 Google Colab Implementation

Because the data is too large, the free Colab doesn't support us to train many times. We call the T100 graphics card to analyze and train in batches when the number of people using the Colab is small, and we also multitask in the same Colab notebook to fully utilize its performance.

The code provided here closely reflects our process, with only minor modifications to ensure reproducible results. Due to our multitasking approach in Colab, some execution paths may vary slightly.

Here is the [Colab Notebook](https://colab.research.google.com/drive/1U6r1HXEAvhQQuxfe2vdppPG_fQSRPjcb?usp=sharing), which is also provided in this repository as `BLAST_Colab.ipynb`.

---

## 🧬 Data Sources

The sequences used in this analysis were obtained from the NCBI GenBank database:
- [Human mitochondrion genome (NC_012920)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920)
- [E. coli genome (NC_000913.3)](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
- [Mouse mitochondrion genome (NC_005089.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_005089.1)
- [Human BRCA1 gene (NG_005905.2)](https://www.ncbi.nlm.nih.gov/nuccore/NG_005905.2)

Due to the large size of these files, they are not included directly in the repository.

---

## 📌 Conclusion

This tool helps explore when BLASTn or BLASTx performs better depending on sequence length and complexity. It can be easily extended to more sequences or to include BLASTp.

---

## ✍️ Author

CS123A Term Project  
Spring 2025

Ray Zhang: rui.zhang@sjsu.edu  
Jisheng Jiang: jisheng.jiang@sjsu.edu