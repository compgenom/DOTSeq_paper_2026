FROM bioconductor/bioconductor_docker:3.22-R-4.5.2

RUN curl -OL https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/miniconda2 && \
    /opt/miniconda2/bin/pip install matplotlib==2.2.5 scipy==1.2.1 statsmodels==0.10.2 && \
    git clone https://github.com/ratschlab/RiboDiff.git /opt/RiboDiff && \
    R -e "BiocManager::install(c('compgenom/DOTSeq', 'devtools', 'optparse', 'testthat', 'BiocCheck', 'BiocStyle', 'rmarkdown', 'knitr', 'DHARMa', 'biomaRt', 'anota2seq', 'eulerr', 'ggplot2', 'ggsignif', 'magick', 'org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg38', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'TxDb.Dmelanogaster.UCSC.dm3.ensGene', 'pasillaBamSubset'), ask = FALSE)" && \
    R -e "options(repos = BiocManager::repositories()); BiocManager::install('iCOBRA', ask = FALSE)" && \
    git clone https://github.com/xryanglab/xtail.git && \
    sed -i 's/\[at\]/@/g' xtail/DESCRIPTION && \
    R -e "devtools::install_local('./xtail')" && \
    git clone https://github.com/smithlabcode/riborex.git && \
    R -e "devtools::install_local('./riborex')" && \
    rm -rf ./xtail ./riborex Miniconda2-latest-Linux-x86_64.sh
