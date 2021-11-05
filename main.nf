#!/usr/bin/env nextflow

/*
 * =======================================================================
 * 
 *		Bidirectional Transcription Detection Pipeline
 * 		
 * =======================================================================
 * 
 * This Source Code is distributed under the GPL-3.0 License
 */


/* 
 * 'Bidirectional-Flow' - A Nextflow pipeline for detecting regions of transcription
 * from Nascent RNA sequencing  
 *
 * The pipeline identifies regions of transcription with FStitch, TFit and/or dREG
 * 
 * =============
 * Authors
 * =============
 * Rutendo F. Sigauke : rutendo.sigauke@cuanschutz.edu
 * Lynn Sanford : lynn.sanford@colorado.edu
 */

def helpMessage() {
    log.info"""
    =========================================
     BidirectionalFlow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile slurm --crams '/project/*.sorted.cram' --workdir '/project/tempfiles' --outdir '/project/'
    Required arguments:
         -profile                      Configuration profile to use. <genome_user>
         --crams                       Directory pattern for cram files: /project/*.sorted.cram (Required if --bams not specified).
         --bams                        Directory pattern for bam files: /project/*.sorted.bam (Required if --crams not specified).
         --workdir                     Nextflow working directory where all intermediate files are saved.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.

    Analysis Options:
        --gene_count                   Run featureCounts to obtain stranded and unstranded gene counts over an annotation.
        --fstitch                      Run FStitch. If used, you must also specify FS_path and FS_train params.
        --tfit                         Run Tfit bidir. If used, you must also specify the Tfit_path parameter.
        --tfit_prelim		       Run Tfit bidir. If used, you must also specify the Tfit_path parameter. Not compatible with --prelim_files flag.
        --tfit_model		       Run Tfit full model. If used, must specify the Tfit path parameter AND have prelim files from --tfit_prelim process or previous run via the --prelim_files flag. Not compatible with --tfit flag.
        --prelim_files		       Directory pattern for tfit prelim files: /project/*-1_prelim_bidir_hits.bed (required for --tfit_model if --tfit_prelim is not also specified)
        --dreg                         Produce bigwigs formatted for input to dREG.
    """.stripIndent()
}

// Configure Variables

params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.tfit_run = "$baseDir/bin/tfit_run.sh"
software_versions = Channel.create()

import java.text.SimpleDateFormat
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd")
output_date =  sdf.format(date)

String output_date = new java.text.SimpleDateFormat("yyMMdd").format(new Date())

// Header log info
log.info """=======================================================
Bidirectional-Flow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'Bidirectional-Flow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = workflow.runName
if(params.crams) summary['Crams']            = params.crams
if(params.bams) summary['Bams']              = params.bams
summary['Genome Ref']       = params.genome
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output dir']       = params.outdir
summary['FStitch']          = params.fstitch ? 'YES' : 'NO'
summary['Tfit']             = params.tfit ? 'YES' : 'NO'
summary['Tfit prelim']      = params.tfit_prelim ? 'YES' : 'NO'
summary['Tfit model']       = params.tfit_model ? 'YES' : 'NO'
summary['Tfit prelim files']	= params.prelim_files ? 'YES' : 'NO'
summary['dREG']             = params.dreg ? 'YES' : 'NO'
summary['Gene counting']    = params.gene_count ? 'YES' : 'NO'
if(params.fstitch)summary['FStitch dir']      = params.fstitch_path
if(params.fstitch)summary['FStitch train']    = params.fstitch_train
if(params.tfit)summary['Tfit dir']      = params.tfit_path
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="


println "\nSTARTING PIPELINE"

// PART 0: Outputting software versions

println "[Log 0]: Getting software versions"

process get_software_versions {
    validExitStatus 0,1,127
    time '1h'

    output:
    stdout into software_versions

    script:
    """
    printf "bidirectionalflow_version: %s\n" ${params.version}
    printf "nextflow_version: %s\n" ${workflow.nextflow.version}
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}')
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}')
    printf "openmpi_version: %s\n" \$(ompi_info | head -2 | tail -1 | awk '{print \$NF}')
    printf "gcc_version: %s\n" \$(gcc --version | head -1 | awk '{print \$NF}')
    printf "fstitch_version: %s\n" \$(${params.fstitch_path} train --version | head -1)
    printf "tfit_version: %s\n" \$(${params.tfit_path} model --version | head -1)
    printf "r_version: %s\n" \$(R --version | head -1 | awk '{print \$3}')
    printf "rsubread_version: %s\n" \$(Rscript -e 'library("Rsubread");packageVersion("Rsubread")' 2>&1 | tail -1 | awk '{print \$NF}')
    printf "boost_version: %s\n" \$(ls -d /Users/\$USER/.local/boost* | head -1 | awk -F "_" '{print \$(NF-2)"."\$(NF-1)"."\$(NF)}')
    printf "dreg_version: %s\n" \$(Rscript -e 'library("dREG");packageVersion("dREG")' 2>&1 | tail -1 | awk '{print \$NF}')
    printf "pipeline_hash: %s\n" ${workflow.scriptId}
    """
}

software_versions.collectFile(name: "software_versions_bidir_${output_date}_${workflow.runName}.yaml", storeDir: "${params.outdir}/pipeline_info")

println "[Log 0]: Software versions complete"


// PART 1: Converting cram files to bam files

if (params.crams) {
  println "[Log 1]: Converting CRAM files to BAM files"
  println "[Log 1]: Genome file being used ..... $params.genome "
  println "[Log 1]: Cram file directory ........ $params.crams"
  println "[Log 1]: Working directory ... $params.workdir"
  println "[Log 1]: Output directory ... $params.outdir"
  cramfiles = Channel
                  .fromPath(params.crams)
		  .map { file -> tuple(file.baseName, file)}		  

  process cram_to_bam {
     cpus 16
     queue 'short'
     memory '5 GB'
     time '1h30m'
     tag "$prefix"

     input:
     set val(prefix),file(cram) from cramfiles

     output:
     set val(prefix), file("${prefix}.sorted.bam") into sorted_bam_file, bam_for_dreg, bam_for_gene_counting, bam_for_bidir_counting
     set val(prefix), file("${prefix}.sorted.bam.bai") into sorted_bam_index, bam_index_for_dreg, bam_index_for_gene_counting, bam_index_for_bidir_counting

     module 'samtools'
     script:
     """
     samtools view -@ 16 -b -1 -T ${params.genome} ${cram} > ${prefix}.sorted.bam
     samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
     """
  }
} else {

  sorted_bam_file = Channel
                  .fromPath(params.bams)
                  .map { file -> tuple(file.baseName, file)}

}

process bam_conversion_tfit {
   cpus 16
   queue 'short'
   memory '5 GB'
   time '1h'
   tag "$prefix"

   when:
   params.tfit || params.fstitch || params.tfit_model || params.tfit_prelim

   input:
   set val(prefix),file(bam) from sorted_bam_file

   output:
   set val(prefix), file("${prefix}.mmfilt.sorted.bam") into bam_for_tfit
   set val(prefix), file("${prefix}.mmfilt.sorted.bam.bai") into bam_index_for_tfit

   module 'samtools'
   script:
   """
   samtools view -@ 16 -h -q 1 ${bam} | \
       grep -P '(NH:i:1|^@)' | \
       samtools view -h -b > ${prefix}.mmfilt.sorted.bam
   samtools index ${prefix}.mmfilt.sorted.bam ${prefix}.mmfilt.sorted.bam.bai
   """
}

println "[Log 1]: Bam files are ready\n"


// PART 2: Generate bedgraphs

process bedgraphs {
    println "[Log 2]: Generating BEDGRAPHS for TFit and FStitch"
    println "[Log 2]: Genome information ..... $params.genome "
    println "[Log 2]: Chromosome Sizes ....... $params.chrom_sizes"

    validExitStatus 0,143
    tag "$prefix"
    memory '40 GB'
    queue 'short'
    time '4h'
    
    publishDir "${params.outdir}/bedgraphs", mode: 'copy', pattern: "${prefix}.bedGraph"

    when:
    params.tfit || params.fstitch || params.tfit_model || params.tfit_prelim

    input:
    set val(prefix), file(bam_file) from bam_for_tfit
    set val(prefix), file(bam_indices) from bam_index_for_tfit

    output:
    set val(prefix), file("${prefix}.pos.bedGraph") into pos_non_normalized_bedgraphs, pos_fstitch_bg
    set val(prefix), file("${prefix}.neg.bedGraph") into neg_non_normalized_bedgraphs, neg_fstitch_bg
    set val(prefix), file("${prefix}.bedGraph") into non_normalized_bedgraphs, fstitch_bg, tfit_bg, prelimtfit_bg, modeltfit_bg, nqc_bg

    script:
    if (params.singleEnd) {
    """
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${bam_file} \
        > ${prefix}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${bam_file} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.neg.bedGraph
    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    """
    } else {
    """
    samtools view \
        -h -b -f 0x0040 \
        ${bam_file} \
        > ${prefix}.first_pair.bam

    samtools view \
        -h -b -f 0x0080 \
        ${bam_file} \
        > ${prefix}.second_pair.bam

    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        > ${prefix}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.first_pair.neg.bedGraph

    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | sortBed \
        > ${prefix}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${prefix}.second_pair.neg.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.pos.bedGraph ${prefix}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.pos.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.neg.bedGraph ${prefix}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.neg.bedGraph

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    """
    }
 }


println "[Log 2]: Bedgraph files have been generated\n"


// PART 3: Running FStitch

process FStitch {
    println "[Log 3]: Running FStitch"
    println "[Log 3]: FStitch training file .. $params.fstitch_train"
    println "[Log 3]: FStich source code ..... $params.fstitch_path"

    tag "$prefix"
    memory '50 GB'
    queue 'short'
    time '4h'
    validExitStatus 0,134
    errorStrategy 'ignore'

    publishDir "${params.outdir}/fstitch/", mode: 'copy', pattern: "*.hmminfo"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "${prefix}.fstitch_seg.ON.merged.bed"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "${prefix}.fstitch_seg.bed"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "*.fstitch_seg.{pos,neg}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "${prefix}.fstitch_bidir.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "*fstitch_bidir.{short,long}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/hist/", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/fstitch/bidirs/stats/", mode: 'copy', pattern: "*.txt"

    when:
    params.fstitch

    input:
    set val(prefix), file(bg) from fstitch_bg
    set val(prefix), file(pos_bg) from pos_fstitch_bg
    set val(prefix), file(neg_bg) from neg_fstitch_bg

    output:
    file ("*.hmminfo") into fs_train_out
    set val(prefix), file ("*.fstitch_seg.ON.merged.bed") into fs_seg_out
    set val(prefix), file ("*.fstitch_bidir.bed") into fs_bidir_out
    file ("*fstitch_bidir.{short,long}.bed") into fs_bidir_short_long_out
    file ("*.html") into fs_bidir_plot_out
    file ("*.txt") into fs_bidir_stats_out

    script:
    """
    ${params.fstitch_path} train \
        -s + \
        -b ${bg} \
        -t ${params.fstitch_train} \
        -o ${prefix}.fstitch.hmminfo

    ${params.fstitch_path} segment \
        -s + \
        -b ${pos_bg} \
        -p ${prefix}.fstitch.hmminfo \
        -o ${prefix}.fstitch_seg.pos.bed

    ${params.fstitch_path} segment \
        -s - \
        -b ${neg_bg} \
        -p ${prefix}.fstitch.hmminfo \
        -o ${prefix}.fstitch_seg.neg.bed

    cat ${prefix}.fstitch_seg.pos.bed \
        ${prefix}.fstitch_seg.neg.bed \
        | sortBed > ${prefix}.fstitch_seg.bed
	    
    cat ${prefix}.fstitch_seg.bed | \
    	grep ON | \
	bedtools merge -i stdin \
	> ${prefix}.fstitch_seg.ON.merged.bed
    
    bidir \
        -b ${prefix}.fstitch_seg.bed \
        -g ${params.genome_refseq} \
        -o ${prefix}.fstitch_bidir.bed \
        -p \
        -s
    """
}

println "[Log 3]: Done Running FStitch\n"


// PART 4: Running Tfit

process tfit {
    println "[Log 4]: Running Tfit (full)"

    tag "$prefix"
    memory '70 GB'
    time '72h'
    cpus 32
    queue 'long'
    validExitStatus 0

    publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*-1_bidir_predictions.bed"
    publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    publishDir "${params.outdir}/tfit/prelim", mode: 'copy', pattern: "*-1_prelim_bidir_hits.bed"

    when:
    params.tfit

    input:
    set val(prefix), file(bg) from tfit_bg

    output:
    file ("*-1_prelim_bidir_hits.bed") into tfit_full_prelim_out
    file ("*-1_bidir_predictions.bed") into tfit_full_bed_out
    file ("*.tsv") into tfit_full_model_out
    file ("*.log") into tfit_full_logs_out

    script:
        """
	${params.tfit_run} -t ${params.tfit_path} \
			   -c ${params.tfit_config} \
			   -b ${bg} \
			   -p ${prefix} \
			   -n 32

        """
}

println "[Log 4]: Done Running Tfit\n"

// PART 4b: Running individual Tfit components

if (params.tfit == false) {
    if (params.prelim_files) {
        tfit_prelim_out = Channel
            .fromPath(params.prelim_files)
            .map { file -> tuple(file.baseName, file)}
    } else if ((params.tfit_model && params.prelim_files == false) || params.tfit_prelim) {
        process tfit_prelim {
            println "[Log 4b]: Running Tfit prelim"

            tag "$prefix"
            memory '70 GB'
            time '72h'
            cpus 32
            queue 'long'
            validExitStatus 0

            publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"
            publishDir "${params.outdir}/tfit/prelim", mode: 'copy', pattern: "*-1_prelim_bidir_hits.bed"

            input:
            set val(prefix), file(bg) from prelimtfit_bg

            output:
            file ("*-1_prelim_bidir_hits.bed") into tfit_prelim_out
            file ("*.tsv") into prelimtfit_model_out
            file ("*.log") into prelimtfit_logs_out

            script:
                """
                sh ${params.tfit_run} -t ${params.tfit_path} \
                                   -c ${params.tfit_config} \
                                   -b ${bg} \
                                   -p ${prefix} \
                                   -n 32

                """
            println "[Log 4b]: Done Running Tfit prelim\n"
        }
    }

    if (params.tfit_model) {
        process tfit_model {
            println "[Log 4b]: Running Tfit model"

            tag "$prefix"
            memory '70 GB'
            time '72h'
            cpus 32
            queue 'long'
            validExitStatus 0

            publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*-1_bidir_predictions.bed"
            publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"

            input:
            set val(prefix), file(bg) from modeltfit_bg
            set file(prelim) from tfit_prelim_out

            output:
            file ("*-1_bidir_predictions.bed") into tfit_model_bed_out
            file ("*.tsv") into tfit_model_model_out
            file ("*.log") into tfit_model_logs_out

            script:
                """

#                ${params.tfit_run} -t ${params.tfit_path} \
#                                   -c ${params.tfit_config} \
#                                   -b ${bg} \
#                                   -p ${prefix} \
#                                   -n 32

                """
        }
    }

println "[Log 4b]: Done Running Tfit model\n"

}


// PART 5: Preparing bigwig files for dREG

process dreg_prep {
    println "[Log 5]: Generating bigwig files for dREG"

    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$prefix"
    memory '60 GB'
    queue 'short'

    publishDir "${params.outdir}/bigwig/", mode: 'copy', pattern: "*.bw"

    when:
    params.dreg

    input:
    set val(prefix), file(bam_file) from bam_for_dreg
    set val(prefix), file(bam_index) from bam_index_for_dreg

    output:
    set val(prefix), file("*.bw") into dreg_bigwig
    set val(prefix), file("${prefix}.pos.bw") into pos_dreg_bw
    set val(prefix), file("${prefix}.neg.bw") into neg_dreg_bw

    script:
    if (params.singleEnd) {
    """
    echo "Creating BigWigs suitable as inputs to dREG"

    export CRAM_REFERENCE=${params.genome}

    bamToBed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${prefix}.dreg.bed
    sortBed -i ${prefix}.dreg.bed > ${prefix}.dreg.sort.bed

    echo "positive strand processed to bedGraph"

    bedtools genomecov \
            -bg \
            -i ${prefix}.dreg.sort.bed \
            -g ${params.chrom_sizes} \
            -strand + \
            > ${prefix}.pos.bedGraph

    sortBed \
            -i ${prefix}.pos.bedGraph \
            > ${prefix}.pos.sort.bedGraph

    bedtools genomecov \
            -bg \
            -i ${prefix}.dreg.sort.bed \
            -g ${params.chrom_sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${prefix}.neg.bedGraph

    sortBed \
            -i ${prefix}.neg.bedGraph \
            > ${prefix}.neg.sort.bedGraph

    echo "negative strand processed to bedGraph"

    ${params.bedGraphToBigWig} ${prefix}.pos.sort.bedGraph ${params.chrom_sizes} ${prefix}.pos.bw
    ${params.bedGraphToBigWig} ${prefix}.neg.sort.bedGraph ${params.chrom_sizes} ${prefix}.neg.bw

    echo "bedGraph to bigwig done"
    """
  } else {
    // https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh
     println "[Log 5]: TO DO- add PE bigwig file conversion https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh"
}

}

println "[Log 5]: Bigwig files are ready \n"


// PART 6: Running dREG

process dreg_run {
    println "Log[6]: Running dREG"
    println "Log[6]: N.B. Requires GPUs"

    tag "$prefix"
    memory '50 GB'
    time '48h'
    cpus 4
    queue 'titan'
    clusterOptions '--gres=gpu'
    validExitStatus 0

    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*dREG*"

    when:
    params.dreg

    input:
    set val(prefix), file(pos_bw) from pos_dreg_bw
    set val(prefix), file(neg_bw) from neg_dreg_bw

    output:
    set val(prefix), file ("${prefix}.*") into dREG_out

    module 'samtools/1.10'

    script:
        """
        bash ${params.dreg_path} \
	     ${pos_bw} \
	     ${neg_bw} \
	     $prefix \
	     ${params.dreg_train} \
	     4 1 
        """
}

// PART 7: Counting over genes

process gene_count {
   println "[Log 7]: Running FeatureCounts"

    tag "$prefix"
    memory '8 GB'
    time '3h'
    cpus 8
    queue 'short'
    validExitStatus 0

    publishDir "${params.outdir}/featurecounts_genes/", mode: 'copy', pattern: "*gene_counts.txt"

    when:
    params.gene_count

    input:
    set val(prefix), file(bam_file) from bam_for_gene_counting
    set val(prefix), file(bam_indices) from bam_index_for_gene_counting

    output:
    set val(prefix), file ("*gene_counts.txt") into gene_count_out

    script:
    if (params.singleEnd) {
        paired = 'FALSE'
    } else {
        paired = 'TRUE'
    }
    
    """
    #!/usr/bin/env Rscript

    library("Rsubread")

    gtf_table <- read.table("${params.filtered_refseq}")   

    if (${paired} == 'FALSE') {

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=1,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=1,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    } else {

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=2,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=2,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    }

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE) 

    gtf_table <- read.table("${params.trunc_refseq}")
    
    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)
    
    """
}

println "[Log 7]: Done Running FeatureCounts\n"


// PART 7: Counting over bidirectionals

process bidir_count {
   println "[Log 7]: Running FeatureCounts (bidirectionals)"

    tag "$prefix"
    memory '8 GB'
    time '3h'
    cpus 8
    queue 'short'
    validExitStatus 0

    publishDir "${params.outdir}/featurecounts_bidirs/", mode: 'copy', pattern: "*bidir_counts.txt"

    when:
    params.bidir_count

    input:
    set val(prefix), file(bam_file) from bam_for_bidir_counting
    set val(prefix), file(bam_indices) from bam_index_for_bidir_counting

    output:
    set val(prefix), file ("*bidir_counts.txt") into bidir_count_out

    script:
    if (params.singleEnd) {
        paired = 'FALSE'
    } else {
        paired = 'TRUE'
    }

    """

    #!/usr/bin/env Rscript

    library("Rsubread")

## Need to make gtf files from bidir beds ##
    gtf_table <- read.table("${params.filtered_refseq}")

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
##        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
#    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
#    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.bidir_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    """
}

println "[Log 8]: Done Running FeatureCounts (bidirs)\n"


/*
 * Completion report
 */
workflow.onComplete {

    def report_fields = [:]
    report_fields['version'] = params.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    report_fields['summary']['Pipeline repository Git URL'] = workflow.repository ?: 'Not stored'
    report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId ?: 'See hash'
    report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision ?: 'See hash'
    report_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    report_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    report_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/report_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/report_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.txt" )
    output_tf.withWriter { w -> w << report_txt }

    log.info "[Bidirectional-Flow] Pipeline Complete"

}

