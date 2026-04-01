# DOTSeq enables genome-wide detection of differential ORF usage

This repository contains reproducible scripts for the [manuscript](https://doi.org/10.1101/2025.09.24.678201).

All precomputed datasets are available on [Zenodo](https://doi.org/10.5281/zenodo.19266548).

## CONTENTS

- **Apptainer**: [app/](app/) — Build the container from an Apptainer definition file. 
- **Reference Files**: [ref/](ref/) — Contains scripts to download reference annotations and transcript files, as well as metadata for both bulk and single-cell analyses.
- **Bulk Datasets**: [bulk/](bulk/) — Scripts to analyse bulk datasets from [Ly 2024](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808).
- **Single-cell Datasets**: [sc/](sc/) — Quantification, aggregation, and visualization of single-cell Ribo-seq data from [VanInsberghe 2021](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481).
- **Benchmarking**: [benchmarking/](benchmarking/) — Scripts for generating simulated datasets and benchmarking DOTSeq across conditions.

## WORKFLOW

### 1. Clone the repository

```bash
git clone https://github.com/compgenom/DOTSeq_paper_2026.git
cd DOTSeq_paper_2026
```

### 2. Set up the computational environment

To install Apptainer, please refer to the official [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html) and follow the instructions for your operating system.

Build the Apptainer container to install all required dependencies, then enter the container:

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
  -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz \
  -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  -o ref
```

### NOTES
#### Data & directory setup

All files are organised into a predefined directory structure to ensure compatibility with downstream steps:
Running the preprocessing scripts will automatically create the required directories.
If you use a custom directory structure or plan to [skip the download/ preprocessing steps](#Skipping-download-/-preprocessing), update paths accordingly when running scripts.

#### Skipping download / preprocessing 

If you prefer not to run the full pipeline, you can start from different stages by downloading the precomputed files available on [Zenodo](https://doi.org/10.5281/zenodo.19266548) accordingly:

1. Start from aligned BAM files
   - [Bulk](#Analysing-bulk-Ribo-seq):
     `*_Aligned.sortedByCoord.out.bam` → place in `data/bulk`
   - [Single-cell](#Analysing-single-cell-Ribo-seq):
     `*_Aligned.sortedByCoord.out_CB.bam` → place in `data/sc/alignment`
2. Start from ORF read counts
   - [Bulk](#Analysing-bulk-Ribo-seq):
     `bulk.rds` → place in `data/bulk/quantification`
   - [Single-cell](#Analysing-single-cell-Ribo-seq):
     `*_Aligned.sortedByCoord.out_CB_mat.rds` → place in `data/sc/quantification`

### 5. Bulk Ribo-seq analysis

#### Downloading raw reads

For the **Bulk datasets**, only samples corresponding to **cycloheximide (CHX)-treated cells** are downloaded, as these are the conditions analysed in the manuscript. This represents a subset of the full dataset from the original study.

```bash
bash bulk/sra_download_bulk.sh \
  --bulk_rna data/bulk/rna \
  --bulk_ribo data/bulk/ribo \
  --threads 16
```
#### Preprocessing raw reads

Before analysis, raw reads must be quality controlled, trimmed, and aligned to generate alignment files.

```bash
bash bulk/preprocessing_bulk.sh
```

#### Analysing bulk Ribo-seq

Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate gene-level counts, and generate visualisations. 

By default, the pipeline starts by using the `*_Aligned.sortedByCoord.out.bam` files to extract exonic reads using `getExonicReads()` from [DOTSeq](https://doi.org/doi:10.18129/B9.bioc.DOTSeq), producing `*_Aligned.sortedByCoord.out.exonic.sorted.bam`. These filtered BAM files are then used for read counting via [DOTSeq](https://doi.org/doi:10.18129/B9.bioc.DOTSeq)'s `countReads()` function. Precomputed exonic BAM files are also available on [Zenodo](https://doi.org/10.5281/zenodo.19266548).

```bash
Rscript bulk/bulk_analysis.R \
  -ss 1 \
  -gr ref/gr_orfs.rds \
  -bam data/bulk \
  -mat data/bulk/quantification \
  -o results/bulk
```

Alternatively, to start from precomputed ORF read counts:

```bash
Rscript bulk/bulk_analysis.R \
  -ss 2 \
  -mat data/bulk/quantification \
  -o results/bulk
```
   
### 6. Single-cell Ribo-seq analysis

#### Installation 

The alignment files for single-cell datasets were prepared using a modified Nextflow pipeline from [scRiboSeq_manuscript](https://github.com/mvanins/scRiboSeq_manuscript). 

**Important**
If you have been running the bulk analysis inside the Apptainer container, exit the container first using `exit`. Then install Nextflow using Miniforge3:

```bash
# Install Miniforge3 
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Miniforge3-24.7.1-2-Linux-x86_64.sh
bash Miniforge3-24.7.1-2-Linux-x86_64.sh -b -p $HOME/miniforge3

# Initialise Conda 
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"

# Create a conda environment 
conda create --name nf
conda activate nf

# Install Nextflow and Java
conda install -c bioconda -c conda-forge nextflow openjdk=17 -y # Java version that is compatible with the nextflow version used in this pipeline
```

#### Clone the forked repository

The modified Nextflow pipeline can be found in this forked [repository](https://github.com/compgenom/scRiboSeq_manuscript).

```bash
git clone https://github.com/compgenom/scRiboSeq_manuscript.git sc/scRiboSeq_manuscript
```

#### Download the raw sequencing reads

For the **Single-cell datasets**, only data derived from **hTERT-RPE1 cells** are used, consistent with the manuscript.

```bash
bash sc/sra_download_sc.sh \
  --sc_dir sc/scRiboSeq_manuscript/data_processing \
  --threads 16 \
  --jobs 2
```

#### Generate the alignment files

```bash
bash sc/scRiboSeq_manuscript/data_processing/preprocessing_sc.sh
```

#### Analysing single-cell Ribo-seq

After obtaining the alignment files, re-launch the Apptainer container `apptainer shell app/dotseq.sif` and run the analysis script to reproduce the manuscript figures.
By default, the script starts from the read counting step. The UMI threshold used in the manuscript is 100, but this can be adjusted.

```bash
Rscript sc/sc_analysis.R \
  -ss 1 \
  -gr ref/gr_orfs.rds \
  -bam data/sc/alignment \
  -mat data/sc/quantification \
  -umi 100 \
  -o results/sc
```

Alternatively, to start from precomputed ORF read counts:

```bash
Rscript sc/sc_analysis.R \
  -ss 2 \
  -gr ref/gr_orfs.rds \
  -mat data/sc/quantification \
  -umi 100 \
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

## REFERENCE 
Jimmy Ly, Kehui Xiang, Kuan-Chung Su, Gunter B Sissoko, David P Bartel, Iain M Cheeseman. (2024).
Nuclear release of eIF1 restricts start-codon selection during mitosis.
Nature. DOI: https://doi.org/10.1038/s41586-024-08088-3

Michael VanInsberghe, Jeroen van den Berg, Amanda Andersson-Rolf, Hans Clevers, Alexander van Oudenaarden. (2021).
Single-cell Ribo-seq reveals cell cycle-dependent translational pausing.
Nature. DOI: https://doi.org/10.1038/s41586-021-03887-4