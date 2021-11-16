#!/bin/bash -ue
printf "bidirectionalflow_version: %s
" 0.2
    printf "nextflow_version: %s
" 18.10.1
    printf "samtools_version: %s
" $(samtools --version | head -1 | awk '{print $NF}')
    printf "bedtools_version: %s
" $(bedtools --version | head -1 | awk -F " v" '{print $2}')
    printf "openmpi_version: %s
" $(ompi_info | head -2 | tail -1 | awk '{print $NF}')
    printf "gcc_version: %s
" $(gcc --version | head -1 | awk '{print $NF}')
    printf "fstitch_version: %s
" $(null train --version | head -1)
    printf "tfit_version: %s
" $(null model --version | head -1)
    printf "r_version: %s
" $(R --version | head -1 | awk '{print $3}')
    printf "rsubread_version: %s
" $(Rscript -e 'library("Rsubread");packageVersion("Rsubread")' 2>&1 | tail -1 | awk '{print $NF}')
    printf "boost_version: %s
" $(ls -d /Users/$USER/.local/boost* | head -1 | awk -F "_" '{print $(NF-2)"."$(NF-1)"."$(NF)}')
    printf "dreg_version: %s
" $(Rscript -e 'library("dREG");packageVersion("dREG")' 2>&1 | tail -1 | awk '{print $NF}')
    printf "pipeline_hash: %s
" cc8f8880330f5920b8964f739e290f02
