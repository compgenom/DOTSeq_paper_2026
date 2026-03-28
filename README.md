# DOTSeq enables genome-wide detection of differential ORF usage

Reproducible scripts and resources for the manuscript.

---

## Contents

- **Apptainer**: [app/](app/) — Build the container from an Apptainer definition file  
- **Reference Genome**: [ref/ref.sh](ref/ref.sh) — Contains scripts to download reference annotations and transcript files, as well as metadata for both bulk and single-cell analyses  
- **Bulk Datasets**: [bulk/](bulk/) — Scripts for bulk dataset analysis  
- **Single-cell Datasets**: [sc/](sc/) — Quantification, aggregation, and visualization of single-cell Ribo-seq data  
- **Benchmarking**: [benchmarking/](benchmarking/) — Scripts for generating simulated datasets and benchmarking DOTSeq across conditions  

---

## Workflow / How to Start

### 1. Clone the repository

```bash
git clone https://github.com/compgenom/DOTSeq_paper_2026.git
cd DOTSeq_paper_2026
```

### 2. Set up the computational environment
   
Build the Apptainer image to ensure all dependencies are installed:

```bash
apptainer build app/dotseq.sif app/dotseq.def
```
   
### 3. Download reference genome, annotation files and raw read files

```bash
bash ref/ref.sh
bash bulk/sra_downloads_bulk.sh
bash sc/sra_downloads_sc.sh
```
   
### 4. Analysing the bulk datasets
   
Before the analysis, raw reads need to be quality controlled, trimmed, and aligned:

```bash
bash bulk/preprocessing_bulk.sh
```

Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate to gene-level counts, and generate visualizations.
By default, the script uses alignment files `.bam` generated from preprocessing as input for read counting via the `countReads()` function in `DOTSeq`.

```bash
Rscript bulk/bulk_analysis.R \
  -ss 1
  -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz
  -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz
  -gr ref/gr_orfs.rds \
  -bam data/bulk/alignment \
  -o results/bulk
```

Alternatively, to start from precomputed ORF read counts, download the `d.rds` file from Zenodo and run:

```bash
Rscript bulk/bulk_analysis.R \
  -ss 2
  -gr ref/gr_orfs.rds \
  -mat data/bulk/alignment \
  -o results/bulk
```
   
### 5. Analysing the single-cell datasets

To reproduce the manuscript figures, use precomputed alignment files (.bam) available on Zenodo and place them in:

```bash
data/sc/alignment
```

By default, the script starts from the read counting step. The recommended UMI threshold is 100, but this can be adjusted.

```bash
Rscript sc/sc_analysis.R \
  -ss 1 \
  -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz \
  -gr ref/gr_orfs.rds
  -bam data/sc/alignment \
  -mat data/sc/quantification \
  -umi 100 \
  -o results/sc
```

Alternatively, to start from precomputed ORF read counts, download the .rds file from Zenodo (link), place it in:

```bash
data/sc/quantification
```

and run:

```bash
Rscript sc/sc_analysis.R \
  -ss 2
  -gr ref/gr_orfs.rds \
  -mat data/sc/alignment \
  -o results/sc
```
   
### 6. Benchmarking

To generate benchmarking figures, simulated datasets must first be created:

```bash
bash benchmarking/create_simdata.sh
```

Then run the benchmarking pipeline:

```bash
Rscript benchmarking/run_packages.R
Rscript benchmarking/plotting.R
```