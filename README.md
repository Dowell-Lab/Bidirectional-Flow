# Bidirectional-Flow
Nextflow pipeline for detecting bidirectional transcripts and counting reads from nascent RNA sequencing data.

This pipeline is built on top of and is an expansion of the Nascent-Flow pipeline (https://github.com/Dowell-Lab/Nascent-Flow). Using dREG, FStitch and/or Tfit, this pipeline identifies regions of transcription fron nascent experiments.

# Installation
 ```
 $ git clone https://github.com/Dowell-Lab/Bidirectional-Flow.git
 $ git checkout Tfit_focus
 ```
 
 
# Edits from Main
1. Running Tfit with 3' bedgraphs. 
    * Why: Tfit is more likely to call bidirectionals with significant background noise when using 3' bedgraphs compared to the full read bedgraphs. 
    * What: Tfit (and only Tfit) will be run with 3prime bedgraphs rather than full-read data. These bedgraphs will be saved if savebg parameter is used.
    * Restrictions: Currently, the edits only allow 3' bedgraphs of single-read data. This means the library must be flipped (rcomp field in the metadata should be 1). 
    * How: --tfit_3prime argument (in help description below)
    * Scripts edited: main_hope.nf
2. Consider 800k+ enhancers identified from DBNascent 2023 (https://doi.org/10.1101/2023.12.07.570626) in the preliminary regions for Tfit. AND consider coverage filters on BOTH strands.
    * Why: Tfit can miss calls so including preliminary regions where we know an enhancer exists, regardless of it is being transcribed, can ensure these regions are considered as possible locations. Similar work was done with the TSS previously.
    * What: Mu regions (defaults are 2kb windows (1kb from mu both directions)) are merged to the original Tfit prelim file, replacing overlaps. Overlapping mu regions are merged. TSS regions then are merged, replacing any prelim regions (including mu regions). Then BOTH strands of regions must have at least 5 reads to be included for downstream Tfit analysis. **Warning: Both strands are currently considered for coverage restrictions regardless of if mu regions are included or not.**
        * *Note*: If not using the double strand coverage, many of these regions can remain due to overlapping genes or noise: in my example, I went from 120k to 43k total prelim regions which saved about 70 hours in downstream Tfit run time.
    * Restrictions: You must provide a bed file of the mus you want to consider. I used hg38_master_qc_gc_len_posneg_gene_filt_divconv_uniqueid_MUforTFIT_08-14-23.bed which can be found in /scratch/Shares/dowell/BidirFlow/hoto7260. I currently have not created parameters outside the python script to parameterize this method further. Defaults are to use 2kb (1kb from mu both directions) windows with no buffers. 
    * How: Add the mu file to your config file for it to be considered. 
    * Scripts edited: prelim_filter_hope.py & main_hope.nf
  

# Requirements

Check each tool for configuration requirements.

- Nextflow (version >= 19.10.0, but the latest versions also appear to cause difficulties)
  You can easily run an older version of nextflow without reinstalling by specifying in the nextflow command, i.e. `NXF_VER=19.10.0 nextflow run ...`

- SAMtools

- BEDtools

- Python3 

- FStitch

  (optional - if running FStitch)

  Install from https://github.com/Dowell-Lab/FStitch

- Tfit

  (optional - if running Tfit)

  Install from https://github.com/Dowell-Lab/Tfit

  Note: FStitch and Tfit also require MPI (Open MPI or MPICH) and GCC for the configuration step. On Fiji Tfit is compiled with Open MPI.

- Rsubread R package

  (optional - if running featureCounts read counting)

  Install as user from Bioconductor https://bioconductor.org/packages/release/bioc/html/Rsubread.html

- dREG

  (optional - if running dREG)

  Install dREG from https://github.com/Danko-Lab/dREG

  Note: Instructions to install the dependencies for dREG can be found on their repository.

  Since rphast is nolonger on CRAN, it can be installed from the tar ball (install.packages("rphast_1.6.11.tar.gz", repos=NULL, type="source"))

  Archived from CRAN https://cran.r-project.org/src/contrib/Archive/rphast/

  Github repo https://github.com/CshlSiepelLab/RPHAST

       - dREG models:

       ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models

       mammals and drosophila https://github.com/Danko-Lab/dREG-Model

## Configuration Files

Each run of Bidirectional-Flow requires a configuration file which specifies genome-specific and user specific variables.

Edit the configuration file to reflect the correct cluster paths (for genome specific variables and cluster-wide executables like Tfit and FStitch on Fiji), and user-specific paths (for dREG installation).

User account-installed libraries such as boost and R packages must be added to the path before running Bidirectional-Flow, using commands such as the following:
  ```
  export R_LIBS_USER=/Users/username/.R/3.6.0/
  export PATH=/Users/username/.local/boost_1_75_0/:"$PATH"
  ```

## FStitch Requirements

FStitch can now optionally be run to segment nascent data into active and inactive regions of transcription and annotate bidirectionals (see https://github.com/Dowell-Lab/FStitch). To run FStitch, you must specify additional parameters in your config file including `FS_path` and `FS_train` which are the full path to the FStitch executable (once compiled) and the training file, respectively. See `example.config` for example parameterization. This option can be executed in the pipeline through the `--fstitch` argument. Please note that the FStitch `bidir` module is in Python3 and must also be pip installed (see Package Requirements).
  ```
  pip3 install fstitch-annotate --user
  ```

## Mapped file inputs

In order to run Tfit, our best practices filter CRAM/BAM files for multimapped reads based on the NH:i: field. The Tfit portion of the pipeline will not work properly for mappers that do not output this field, such as Bowtie2. If you have mapped files that do not contain this tag, you can filter based on the XS:i: tag. This involves switching which version of the script is commented out (using // comment characters) within the `bam_conversion_tfit` process.

## Loading requirements on SLURM

  Below is a summary of all FIJI modules needed to run Bidirectional-Flow. The R version you use depends on which version used for installing packages in user accounts.

  ```
  module load samtools/1.8
  module load bedtools/2.28.0
  module load openmpi/1.6.4
  module load gcc/7.1.0
  module load python/3.6.3
  module load R/3.6.1
  ```

# Running Bidirectional-Flow

## Usage

   The typical command for running the pipeline is as follows:

   ```
   nextflow run /Bidirectional-Flow/main.nf -profile hg38 \
    --crams "/paths/to/*cram" \
    --workdir /processed_data/ \
    --singleEnd \
    --tfit \
    --dreg \
    --savebidirs

   ```

   Arguments (see further details about Tfit and dREG options below):

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
        --savestats                    Saves tfit/dreg bidirectional statistics for all samples (default true)
        --savebam                      Saves sorted bam files (if cram as input)
        --savetfitbam                  Saves multimapped read filtered bamfiles used for tfit
        --savebg                       Saves bedgraph files for tfit and/or dreg (whichever is run)
        --savebw                       Saves dreg bigwig files
        --savebidirs                   Saves bedfiles for promoter/intronic/intergenic bidirectional subsets

    Analysis Options:
        --gene_count                   Run featureCounts to obtain stranded and unstranded gene counts over an annotation.
        --bidir_count                  Run featureCounts to obtain stranded and unstranded read counts over an SAF annotation, such as for Tfit- and/or dREG-derived bidirectionals
        --fstitch                      Run FStitch. If used, you must also specify FS_path and FS_train params.
        --tfit                         Run Tfit bidir and full model. If used, you must also specify the Tfit_path parameter.
        --tfit_3prime                  Use just the 3 prime end bedgraphs for Tfit analysis (should be flipped) **JUST WITH main_hope.nf**
        --tfit_prelim                  Run Tfit bidir. If used, you must also specify the Tfit_path parameter. Not compatible with --prelim_files flag.
        --tfit_model                   Run Tfit full model. If used, must specify the Tfit path parameter AND have prelim files from --tfit_prelim process or previous run via the --prelim_files flag. Not compatible with --tfit flag.
        --tfit_split_model             Run Tfit model with different k values for different size regions (<5kb and 5-10kb)
        --prelim_process               Split regions and filter Tfit prelim files according to current best practices(default true). Must be true if running tfit_split_model
        --prelim_files                 Directory pattern for tfit prelim files: /project/*-1_prelim_bidir_hits.bed (required for --tfit_model if --tfit_prelim is not also specified)
        --dreg                         Produce bigwigs formatted for input to dREG.
        --dreg_results                 Run coverage filtering on existing dREG results

   ```

## Tfit options

Tfit (see https://github.com/Dowell-Lab/Tfit) contains an internal template matching algorithm (the `bidir` module) to generate preliminary regions to model, then will model RNAPII activity using the `model` module.

The two most likely methods of running Tfit are as follows:

   ```
   nextflow run /Bidirectional-Flow/main.nf -profile hg38 \
    --crams "/paths/to/*cram" \
    --workdir /processed_data/ \
    --singleEnd \
    --tfit_split_model
   ```

   ```
   nextflow run /Bidirectional-Flow/main.nf -profile hg38 \
    --crams "/paths/to/*cram" \
    --workdir /processed_data/ \
    --singleEnd \
    --tfit
   ```

Running the `--tfit` argument will run both Tfit modules in sequence. Running the `--tfit_split_model` argument will run both modules, but will run the `model` module in an optimized way (implemented for DBNascent). This is detailed further in the "Split Model Run" section below.

With either flag, in between running the two modules, the pipeline will by default run a prelim region processing script, detailed further in the "Prelim Processing" section below.

The modules may also be run separately with `--tfit_prelim`, which runs the `bidir` module, or `--tfit_model` in conjunction with `--prelim_files`, which runs the `model` module on already-generated prelim files. For the second option, prelim processing is again run by default.

### Prelim Processing

`--prelim_process` is true by default. This informs the pipeline to process prelim files prior to modeling by performing three steps:
- Add 2kb regions around annotated transcription start sites
- Cut large prelim regions into equal-size regions less than 10kb.
- Filter out regions that have less than 5 unique reads, as determined with `bedtools coverage`

Before final model output, all modeled output regions are also filtered for coverage using the same metric. This outputs an additional coverage filtered bedfile.

### Split Model Run

`--tfit_split_model` runs the Tfit `model` module in an optimized way. It must be used with `--prelim_process=true` (default behavior).

This flag runs two separate instances of the Tfit `model` module, specifying the model to look for up to two bidirectionals (k=2) in prelim regions <5kb in length and to look for up to five bidirectionals (k=5) in prelim regions from 5-10kb in length. Without this option (using either the `--tfit` or `--tfit_model`) flags, a single instance of the model module is called with a blanket k=5.

Please be aware that the Tfit model module may take a long time to run (5-70 hrs on mammalian datasets, depending on dataset complexity), regardless of which flag is used. The split modeling process may take less walltime, but may use as many resources due to two parallel modeling processes.

## dREG options

dREG (see https://github.com/Danko-Lab/dREG) is run according to specifications from the Danko lab.

By default, dREG outputs undergo a coverage filtering step identical to that performed on Tfit outputs (see above). This coverage filtering process can be run on existing dREG results using the `--dreg_results` flag to specify paths for output files.

## Credits

* Rutendo Sigauke (): Nextflow base code and pipeline implementation
* Lynn Sanford ([@lynn-sanford](https://github.com/lynn-sanford)): Nextflow pipeline implementation and optimization
* Hope Townsend (hope.townsend@colorado.edu) made changes explained in "Edits from Main" section
