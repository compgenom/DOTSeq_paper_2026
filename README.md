# DOTSeq_paper_2026

Reproducible scripts and resources for the **DOTSeq 2026** manuscript.

---

## Contents

- **Apptainer**: [app/](app/) — Build container from an Apptainer definition file. 
- **Reference Genome**: [ref.sh](ref.sh) — Scripts to download reference genome annotations and transcript fasta files.
- **Bulk Datasets**: [bulk/](bulk/) — Contains all the scripts for bulk datasets.
- **Single-cell Datasets**: [sc/](sc/) — Analysis of single-cell Ribo-seq datasets: quantification, aggregation, visualization.
- **Benchmarking**: [benchmarking/](benchmarking/) — Scripts to generate simulated sequencing datasets and for benchmarking DOTSeq across datasets and conditions.

---

## Workflow / How to Start

1. **Set up the computational environment**
   
   Build the Apptainer image to ensure all dependencies are installed:

   ```bash
   apptainer build app/dotseq.sif app/dotseq.def
   ```
2. **Download reference genome, annotations and SRA files**

   Run the preprocessing scripts for reference data:

   ```bash
   bash ref.sh
   bash bulk/sra_downloads_bulk.sh
   bash sc/sra_downloads_sc.sh
   ```
3. **Preprocess raw sequencing data**
   
   Before any analysis, raw reads need to be quality controlled, trimmed, and aligned:
   
   ```bash
   bash sc/preprocessing_sc.sh # single cell dataset (using Nextflow)
   bash bulk/preprocessing_bulk.sh # bulk dataset
   ```
4. **Run analysis scripts**

   Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate to gene-level counts, and generate visualizations:

   ```bash
   Rscript sc/sc_analysis.R # single cell dataset
   Rscript bulk/bulk_analysis.R # bulk dataset
   ```
5. **Benchmarking**

   Simulated datasets must be generated prior to benchmarking analysis.
   ```bash
   bash benchmarking/create_simdata.sh
   Rscript benchmarking/run_packages.R
   Rscript benchmarking/plotting.R
   ```