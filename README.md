# Bidirectional-Flow
Nextflow pipeline for detecting bidirectional transcripts and counting reads from nascent RNA sequencing data.

This pipeline is built on top of and is an expansion of the Nascent-Flow pipeline (https://github.com/Dowell-Lab/Nascent-Flow). Using dREG, FStitch and/or Tfit, this pipeline identifies regions of transcription fron nascent experiments.

# Installation

  $ git clone https://github.com/Dowell-Lab/Bidirectional-Flow.git 

# Requirements

Check each tool for configuration requirements.

- Nextflow

- SAMtools

- BEDtools

- Python3 

- dREG
Install dREG from https://github.com/Danko-Lab/dREG

Note: Instructions to install the dependencies for dREG can be found on their repository.
Since rphast is nolonger on CRAN, it can be nstalled from the tar ball (install.packages("rphast_1.6.11.tar.gz", repos=NULL, type="source"))

Archived from CRAN https://cran.r-project.org/src/contrib/Archive/rphast/

Github repo https://github.com/CshlSiepelLab/RPHAST

-- dREG models:

ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models

mammals and drosophila https://github.com/Danko-Lab/dREG-Model

- FStitch
Install from https://github.com/Dowell-Lab/FStitch

- Tfit
Install from https://github.com/Dowell-Lab/Tfit

Note: FStitch and Tfit also require MPI (Open MPI or MPICH) and GCC for the configuration step.

# Running Bidirectional-Flow



## Loading requirements on SLURM

  Below is a summary of all FIJI modules needed to run TFEA.
  
  ```
  module load samtools/1.8	
  module load bedtools/2.28.0
  module load openmpi/1.6.4
  module load gcc/7.1.0
  module load python/3.6.3
  ```

## Example run

   ```
   nextflow run /Bidirectional-Flow/main.nf -profile hg38 \
    --crams "/paths/to/*cram" \
    --workdir /processed_data/ \
    --singleEnd \
    --tfit_prelim \
    --tfit_model \
    --dreg

   ```
