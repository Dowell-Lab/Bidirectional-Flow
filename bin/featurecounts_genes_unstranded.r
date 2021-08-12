library("Rsubread")

# Read in arguments
inbam <- commandArgs(trailingOnly = TRUE)[1]
outdir <- commandArgs(trailingOnly = TRUE)[2]
prefix <- commandArgs(trailingOnly = TRUE)[3]
genomegtf <- commandArgs(trailingOnly = TRUE)[4]
if (commandArgs(trailingOnly = TRUE)[5] == 'PE') {
    paired <- TRUE
} else {
    paired <- FALSE
}

# Read counting using featureCounts
fc <- featureCounts(files=inbam,
                    annot.ext=genomegtf,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="gene_length",
                    useMetaFeatures=FALSE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=FALSE,
                    isPairedEnd=paired,
                    strandSpecific=0,
                    nthreads=8)

# Write results -- keep our accession number (GeneID) and length
write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            file=paste0(outdir,"/",prefix,".unstranded.gene_counts.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)

