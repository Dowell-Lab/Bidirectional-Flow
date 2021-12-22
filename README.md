# Bidirectional-Flow
Nextflow pipeline for detecting bidirectional transcripts and counting reads from nascent RNA sequencing data.

This pipeline is built on top of and is an expansion of the Nascent-Flow pipeline (https://github.com/Dowell-Lab/Nascent-Flow). Using dREG, FStitch and/or Tfit, this pipeline identifies regions of transcription fron nascent experiments.

# Installation

 `$ git clone https://github.com/Dowell-Lab/Bidirectional-Flow.git`

# Requirements

Check each tool for configuration requirements.

- Nextflow (version >= 19.10.0)

- SAMtools

- BEDtools

- Python3 

- dREG
  Install dREG from https://github.com/Danko-Lab/dREG

  Note: Instructions to install the dependencies for dREG can be found on their repository.
  Since rphast is nolonger on CRAN, it can be nstalled from the tar ball (install.packages("rphast_1.6.11.tar.gz", repos=NULL, type="source"))

  Archived from CRAN https://cran.r-project.org/src/contrib/Archive/rphast/

  Github repo https://github.com/CshlSiepelLab/RPHAST

       - dREG models:

       ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models

       mammals and drosophila https://github.com/Danko-Lab/dREG-Model

- FStitch
  Install from https://github.com/Dowell-Lab/FStitch

- Tfit
  Install from https://github.com/Dowell-Lab/Tfit

  Note: FStitch and Tfit also require MPI (Open MPI or MPICH) and GCC for the configuration step. On Fiji Tfit is compiled with Open MPI.

# Running Bidirectional-Flow

## Usage

   The typical command for running the pipeline is as follows:
    
   `$ nextflow run main.nf -profile slurm --crams '/project/*.sorted.cram' --workdir '/project/tempfiles' --outdir '/project/'`

   Below are the arguments :

   ```
    Required arguments:
         -profile                      Configuration profile to use. <genome_user>
         --crams                       Directory pattern for cram files: /project/*.sorted.cram (Required if --bams not specified).
         --bams                        Directory pattern for bam files: /project/*.sorted.bam (Required if --crams not specified).
         --workdir                     Nextflow working directory where all intermediate files are saved.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --r1_five_prime                If input file is paired, specifies if read 1 has the 5 prime end (default R2 is five prime, must be manually determined)

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.

    Analysis Options:
        --gene_count                   Run featureCounts to obtain stranded and unstranded gene counts over an annotation.
        --fstitch                      Run FStitch. If used, you must also specify FS_path and FS_train params.
        --tfit                         Run Tfit bidir and full model. If used, you must also specify the Tfit_path parameter.
        --tfit_prelim                  Run Tfit bidir. If used, you must also specify the Tfit_path parameter. Not compatible with --prelim_files flag.
        --tfit_model                   Run Tfit full model. If used, must specify the Tfit path parameter AND have prelim files from --tfit_prelim process or previous run via the --prelim_files flag. Not compatible with --tfit flag.
        --tfit_split_model             Run Tfit model with different k values for different size regions (<5kb and 5-10kb)
        --prelim_process               Split regions and filter Tfit prelim files according to current best practices. Automatically run with tfit_split_model
        --prelim_files                 Directory pattern for tfit prelim files: /project/*-1_prelim_bidir_hits.bed (required for --tfit_model if --tfit_prelim is not also specified)
        --dreg                         Produce bigwigs formatted for input to dREG.

   ```

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
    --tfit \
    --prelim_process \
    --dreg

   ```
