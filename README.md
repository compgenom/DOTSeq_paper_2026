# DOTSeq enables genome-wide detection of differential ORF usage

Reproducible scripts and resources for the manuscript.

---

## Contents

- **Apptainer**: [app/](app/) — Build container from an Apptainer definition file. 
- **Reference Genome**: [ref/ref.sh](ref/ref.sh) — A folder containing the scripts to download the reference annotation and transcript files. It also contains the metadata for both bulk and single-cell analysis. 
- **Bulk Datasets**: [bulk/](bulk/) — Contains all the scripts for bulk datasets.
- **Single-cell Datasets**: [sc/](sc/) — Analysis of single-cell Ribo-seq datasets: quantification, aggregation, visualization.
- **Benchmarking**: [benchmarking/](benchmarking/) — Scripts to generate simulated sequencing datasets and for benchmarking DOTSeq across datasets and conditions.

---

## Workflow / How to Start

1. **Clone the repository**

   ```bash
   git clone https://github.com/compgenom/DOTSeq_paper_2026.git
   cd DOTSeq_paper_2026
   ```
2. **Set up the computational environment**
   
   Build the Apptainer image to ensure all dependencies are installed:

   ```bash
   apptainer build app/dotseq.sif app/dotseq.def
   ```
   
3. **Download reference genome, annotation files and raw read files**

   ```bash
   bash ref/ref.sh
   bash bulk/sra_downloads_bulk.sh
   bash sc/sra_downloads_sc.sh
   ```
   
4. **Preprocess raw sequencing data**
   
   Before any analysis, raw reads need to be quality controlled, trimmed, and aligned:
   
   ```bash
   bash bulk/preprocessing_bulk.sh # bulk dataset
   ```

   The full preprocessing steps for the single-cell raw reads will be provided in a seperate repository.
   
5. **Bulk analysis**

   Once preprocessing is complete, run the analysis scripts to quantify ORFs, aggregate to gene-level counts, and generate visualizations.
   By default, the script will take the alignment files (.bam) generated form the preprocessing step as input to do read counting using the countReads() function in DOTSeq:
   
   ```bash
   Rscript bulk/bulk_analysis.R \
      -ss 1
      -a ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz
      -s ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz
      -gr ref/gr_orfs.rds \
      -bam data/bulk/alignment \
      -o results/bulk
   ```

   Alternatively, if you wish to start from the pre-computed ORF read counts, you may download the d.rds file from Zenodo. (link) and run the following command:
   
   ```bash
   Rscript bulk/bulk_analysis.R \
      -ss 2
      -gr ref/gr_orfs.rds \
      -mat data/bulk/alignment \
      -o results/bulk
   ```
6. **Single-cell dataset analysis**
   
   To reproduce the manuscript figures, use the precomputed alignment files (.bam) available on Zenodo and place them in:

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

   Alternatively, if you wish to start from the pre-computed ORF read counts, you may download the .rds file from Zenodo (link), place them in:

   ```bash
   data/sc/quantification
   ```

   and run the following command:
   
   ```bash
   Rscript sc/sc_analysis.R \
      -ss 2
      -gr ref/gr_orfs.rds \
      -mat data/sc/alignment \
      -o results/sc
   ```
   
7. **Benchmarking**

   To generate the figures for the benchmarking analysis, simulated datasets must be generated prior to benchmarking analysis with:
   
   ```bash
   bash benchmarking/create_simdata.sh
   ```

   and then followed by the main pipeline to benchmark DOTSeq against other tools.
   
   ```bash
   Rscript benchmarking/run_packages.R
   Rscript benchmarking/plotting.R
   ```