# DOTSeq_paper_2026

Reproducible scripts and resources for the **DOTSeq 2026** manuscript.

---

## Contents

- **Docker**: [docker/](docker/) — Dockerfile to build the DOTSeq computational environment.
- **Reference Genome**: [scripts/preprocessing/ref.sh](scripts/preprocessing/ref.sh) — Scripts to download reference genome annotations and transcript fasta files.
- **Preprocessing**: [scripts/preprocessing/](scripts/preprocessing/) — Upstream read processing: QC, adapter trimming, alignment.
- **Analysis**: [scripts/analysis/](scripts/analysis/) — Analysis of single-cell Ribo-seq datasets: quantification, aggregation, visualization.
- **Simulation**: [scripts/simulation/](scripts/simulation/) — Scripts to generate simulated sequencing datasets.
- **Benchmarking**: [scripts/benchmarking/](scripts/benchmarking/) — Scripts for benchmarking DOTSeq across datasets and conditions.

---

## Workflow / How to Start

1. **Set up the computational environment**  
   Build the Docker image to ensure all dependencies are installed:

   ```bash
   docker build -t dotseq:2026 docker/
   ```
2. **Download reference genome and annotations**
    Run the preprocessing scripts for reference data:

   ```bash
   bash scripts/preprocessing/ref.sh
   ```
3. **Preprocess raw sequencing data**
    Before any analysis, raw reads need to be quality controlled, trimmed, and aligned:
   
   ```bash
   bash scripts/preprocessing/preprocessing_sc.sh # single cell dataset (using Nextflow)
   bash scripts/preprocessing/preprocessing_bulk.sh # bulk dataset
   ```
4. **Run analysis scripts**
   Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate to gene-level counts, and generate visualizations:

   ```bash
   Rscript scripts/analysis/sc_analysis.R # single cell dataset
   Rscript scripts/analysis/bulk_analysis.R # bulk dataset
   ```
5. **Benchmarking**
   Simulated data is used in the paper for benchmarking analysis. This is the script used to generate simulated datasets:
   ```bash
   bash scripts/simulation/create_simdata.sh
   ```
   Benchmark DOTSeq performance using simulated datasets:
   ```bash
   Rscript scripts/benchmarking/run_packages.R
   Rscript scripts/benchmarking/plotting.R
   ```