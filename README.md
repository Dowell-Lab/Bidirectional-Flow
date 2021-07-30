# Bidirectional-Flow
Nextflow pipeline for detecting bidirectional transcripts from nascent RNA sequencing data.

This pipeline is built on top of and is an expansion of the Nascent-Flow pipeline (https://github.com/Dowell-Lab/Nascent-Flow). Using dREG, FStitch and/or Tfit, this pipeline identifies regions of transcription fron nascent experiments.

Requirements
dREG
Install dREG from https://github.com/Danko-Lab/dREG

Note: Instructions to install the dependencies for dREG can be found on the repository
Since rphast is nolonger on CRAN, it can be nstalled from the tar ball (install.packages("rphast_1.6.11.tar.gz", repos=NULL, type="source"))
Archived from CRAN https://cran.r-project.org/src/contrib/Archive/rphast/
Github repo https://github.com/CshlSiepelLab/RPHAST
dREG models

ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models
mammals and drosophila https://github.com/Danko-Lab/dREG-Model
FSticth
Install from https://github.com/Dowell-Lab/FStitch
TFit
Install from https://github.com/Dowell-Lab/Tfit
