# DOTSeq enables genome-wide detection of differential ORF usage

This repository contains reproducible scripts for the [manuscript](https://doi.org/10.1101/2025.09.24.678201).

All precomputed datasets are available on [Zenodo](https://doi.org/10.5281/zenodo.19266548).

## CONTENTS

- **Apptainer**: [app/](app/) — Build the container from an Apptainer definition file  
- **Reference Genome**: [ref/](ref/) — Contains scripts to download reference annotations and transcript files, as well as metadata for both bulk and single-cell analyses  
- **Bulk Datasets**: [bulk/](bulk/) — Scripts to analyse bulk datasets from [Ly 2024](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
- **Single-cell Datasets**: [sc/](sc/) — Quantification, aggregation, and visualization of single-cell Ribo-seq data from [VanInsberghe 2021](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481) 
- **Benchmarking**: [benchmarking/](benchmarking/) — Scripts for generating simulated datasets and benchmarking DOTSeq across conditions  

---

## WORKFLOW

### 1. Clone the repository

```bash
git clone https://github.com/compgenom/DOTSeq_paper_2026.git
cd DOTSeq_paper_2026
```

### 2. Set up the computational environment
   
Build the Apptainer image to ensure all dependencies are installed and enter the apptainer image:

```bash
apptainer build app/dotseq.sif app/dotseq.def
apptainer shell app/dotseq.sif
```

### 3. Download reference genome, annotation files and raw read files

```bash
bash ref/ref.sh
bash bulk/sra_downloads_bulk.sh
bash sc/sra_downloads_sc.sh
```
   
### 4. Analysing the bulk datasets
   
Before the analysis, raw reads need to be quality controlled, trimmed, and aligned to produce an alignment file.

```bash
bash bulk/preprocessing_bulk.sh
```
These `.bam` alignment files are also available on [Zenodo](https://doi.org/10.5281/zenodo.19266548), if you wish to skip the preprocessing step.

Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate to gene-level counts, and generate visualizations.
By default, the script uses alignment files `.bam` generated from preprocessing as input for read counting via the `countReads()` function in [DOTSeq](https://github.com/compgenom/DOTSeq/tree/main).

```bash
Rscript bulk/bulk_analysis.R \
  -ss 1
  -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz
  -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz
  -gr ref/gr_orfs.rds \
  -bam data/bulk/alignment \
  -mat data/bulk/quantification \
  -o results/bulk
```

Alternatively, to start from precomputed ORF read counts, download the `bulk.rds` file from [Zenodo](https://doi.org/10.5281/zenodo.19266548), place it in:

```bash
data/bulk/quantification
```

and run:

```bash
Rscript bulk/bulk_analysis.R \
  -ss 2
  -gr ref/gr_orfs.rds \
  -mat data/bulk/quantification/bulk.rds \
  -o results/bulk
```
   
### 5. Analysing the single-cell datasets

To reproduce the manuscript figures, use precomputed alignment files `.bam` available on [Zenodo](https://doi.org/10.5281/zenodo.19266548) and place them in:

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

Alternatively, to start from precomputed ORF read counts, download the `.rds` file from [Zenodo](https://doi.org/10.5281/zenodo.19266548), place it in:

```bash
data/sc/quantification
```

and run:

```bash
Rscript sc/sc_analysis.R \
  -ss 2
  -gr ref/gr_orfs.rds \
  -mat data/sc/quantification \
  -umi 100 \
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

## CONTACTS AND BUG REPORTS
- Chun Shen Lim: 
chunshen [dot] lim [at] otago [dot] ac [dot] nz
- Gabrielle Chieng: 
gabrielle [dot] chieng [at] postgrad [dot] otago [dot] ac [dot] nz

## CITATION
compgenom. DOTSeq: Genome-wide detection of differential ORF usage. 2026. GitHub, https://github.com/compgenom/DOTSeq.