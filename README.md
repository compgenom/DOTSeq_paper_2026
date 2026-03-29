# DOTSeq enables genome-wide detection of differential ORF usage

This repository contains reproducible scripts for the [manuscript](https://doi.org/10.1101/2025.09.24.678201).

All precomputed datasets are available on [Zenodo](https://doi.org/10.5281/zenodo.19266548).

## CONTENTS

- **Apptainer**: [app/](app/) — Build the container from an Apptainer definition file  
- **Reference sequence and annotation**: [ref/](ref/) — Contains scripts to download reference annotations and transcript files and generate ORF-level annotation, as well as metadata for both bulk and single-cell analyses
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

### 3. Download reference sequence and annotation files
This step retrieves all reference files required for both bulk and single-cell analyses. 
The following scripts will download the required reference sequence and annotation files into the appropriate directories:

```bash
bash ref/ref.sh
```

### 4. Generate the ORF-level annotation using DOTSeq's `getORFs()` function

```bash
Rscript ref/orf_annotation.R \
  -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz \
  -o ref
```

### 5.Bulk datasets

#### Downloading raw reads

For the **Bulk datasets**, only samples corresponding to **cycloheximide (CHX)-treated cells** are downloaded, as these are the conditions analysed in the manuscript. This represents a subset of the full dataset from the original study.

```bash
bash bulk/sra_download.sh \
  --bulk_rna data/bulk/rna \
  --bulk_ribo data/bulk/ribo \
  --threads 16 \
```
#### Preprocessing raw reads

Before analysis, raw reads must be quality controlled, trimmed, and aligned to generate alignment files. Precomputed alignment files (`Aligned.sortedByCoord.out.bam`) are available on [Zenodo](https://doi.org/10.5281/zenodo.19266548) if you wish to skip this preprocessing step.

```bash
bash bulk/preprocessing_bulk.sh
```

#### Analysing bulk datasets

Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate gene-level counts, and generate visualisations.
By default, the pipeline uses the `Aligned.sortedByCoord.out.bam` files to extract exonic reads using `getExonicReads()` from [DOTSeq](https://github.com/compgenom/DOTSeq/tree/main), producing `_Aligned.sortedByCoord.out.exonic.sorted.bam`. These filtered BAM files are then used for read counting via [DOTSeq](https://github.com/compgenom/DOTSeq/tree/main)'s `countReads()` function. Precomputed exonic BAM files are also available on [Zenodo](https://doi.org/10.5281/zenodo.19266548).

```bash
Rscript bulk/bulk_analysis.R \
  -ss 1 \
  -gr ref/gr_orfs.rds \
  -bam data/bulk/alignment \
  -mat data/bulk/quantification \
  -o results/bulk
```

Alternatively, to start from precomputed ORF read counts, download the `bulk.rds` file from [Zenodo](https://doi.org/10.5281/zenodo.19266548), place it in `data/bulk/quantification` and run:

```bash
Rscript bulk/bulk_analysis.R \
  -ss 2 \
  -gr ref/gr_orfs.rds \
  -mat data/bulk/quantification \
  -o results/bulk
```

##### Notes on directory settings
All files will be organised into the predefined directory structure (e.g., `ref/`, and `data/bulk/`) to ensure compatibility with downstream preprocessing and analysis steps. 
If you run the preprocessing steps, the required directory structure (e.g., `data/bulk/`) will be created automatically.
If you prefer to use a different directory structure, or if you choose to skip the preprocessing step and download precomputed files from [Zenodo](https://doi.org/10.5281/zenodo.19266548), you will need to create the appropriate directories manually before placing the files. Please also modify the paths when running the scripts accordingly.

   
### 6. Single-cell datasets

For the **Single-cell datasets**, only data derived from **hTERT-RPE1 cells** are used for downstream analysis, consistent with the conditions presented in the manuscript.

The preprocessing step to generate the alignment files for single-cell datasets was done using a modified Nextflow pipeline from [scRiboSeq_manuscript](https://github.com/mvanins/scRiboSeq_manuscript). This modified pipeline will be provided in a separate repository.

To reproduce the manuscript figures, download the precomputed alignment files `_Aligned.sortedByCoord.out_CB.bam` from on [Zenodo](https://doi.org/10.5281/zenodo.19266548) and place them in `data/sc/alignment`.

By default, the script starts from the read counting step. The recommended UMI threshold is 50, but this can be adjusted.

```bash
Rscript sc/sc_analysis.R \
  -ss 1 \
  -gr ref/gr_orfs.rds \
  -bam data/sc/alignment \
  -mat data/sc/quantification \
  -umi 50 \
  -o results/sc
```

Alternatively, to start from precomputed ORF read counts, download the `_Aligned.sortedByCoord.out_CB_mat.rds` and `gr_orfs.rds` file from [Zenodo](https://doi.org/10.5281/zenodo.19266548), place them in `data/sc/quantification` and `ref/` respectively, and run:

```bash
Rscript sc/sc_analysis.R \
  -ss 2 \
  -gr ref/gr_orfs.rds \
  -mat data/sc/quantification \
  -umi 50 \
  -o results/sc
```
   
### 7. Benchmarking

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
Chun Shen Lim, Gabrielle S. W. Chieng. (2026). 
DOTSeq enables genome-wide detection of differential ORF usage. 
BioRxiv. DOI: https://doi.org/10.1101/2025.09.24.678201