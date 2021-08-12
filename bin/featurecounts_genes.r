library("Rsubread")

# Read in arguments
inbam <- commandArgs(trailingOnly = TRUE)[1]
outdir <- commandArgs(trailingOnly = TRUE)[2]
prefix <- commandArgs(trailingOnly = TRUE)[3]
genomegtf <- commandArgs(trailingOnly = TRUE)[4]

# Read counting using featureCounts
fc <- featureCounts(files=inbam,
                    annot.ext=genomegtf,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="whole_gene",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=FALSE,
                    isPairedEnd=FALSE,
                    nthreads=8)

# Write results -- keep our accession number (GeneID) and length
write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            paste0(outdir,prefix,file=".gene_counts.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)

