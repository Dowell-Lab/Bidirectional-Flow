#!/usr/bin/env nextflow

/*
 * =======================================================================
 * 
 *		Nascent Transcription Detection Pipeline
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
 * Rutendo F. Sigauke
 * rutendo.sigauke@cuanschutz.edu
 */


// Configure Variables

params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"

println "\nSTARTING PIPELINE"

// PART 1: Converting cram files to bam files

if (params.crams) {
  println "[Log 1]: Converting CRAM files to BAM files"
  println "[Log 1]: Genome file being used ..... $params.genome "
  println "[Log 1]: Cram file directory ........ $params.crams"
  println "[Log 1]: Output/Working directory ... $params.workdir"
  cramfiles = Channel
                  .fromPath(params.crams)
		  .map { file -> tuple(file.baseName, file)}
		  

  process cram_to_bam {
     cpus 16
     queue 'short'
     memory '50 GB'
     time '1h30m'
     tag "$prefix"

     input:
     set val(prefix),file(cram) from cramfiles

     output:
     set val(prefix), file("${prefix}.sorted.bam") into sorted_bam_file, bam_for_dreg
     set val(prefix), file("${prefix}.sorted.bam.bai") into sorted_bam_index, bam_index_for_dreg

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
    publishDir "${params.workdir}/bedgraphs", mode: 'copy', pattern: "${prefix}.bedGraph"

    input:
    set val(prefix), file(bam_file) from sorted_bam_file
    set val(prefix), file(bam_indices) from sorted_bam_index

    output:
    set val(prefix), file("${prefix}.pos.bedGraph") into pos_non_normalized_bedgraphs, pos_fstitch_bg
    set val(prefix), file("${prefix}.neg.bedGraph") into neg_non_normalized_bedgraphs, neg_fstitch_bg
    set val(prefix), file("${prefix}.bedGraph") into non_normalized_bedgraphs, fstitch_bg, tfit_bg, prelimtfit_bg, nqc_bg

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

    publishDir "${params.workdir}/fstitch/", mode: 'copy', pattern: "*.hmminfo"
    publishDir "${params.workdir}/fstitch/segment/", mode: 'copy', pattern: "*.fstitch_seg.bed"
    publishDir "${params.workdir}/fstitch/bidirs/", mode: 'copy', pattern: "${prefix}.fstitch_bidir.bed"
    publishDir "${params.workdir}/fstitch/bidirs/", mode: 'copy', pattern: "*fstitch_bidir.{short,long}.bed"
    publishDir "${params.workdir}/fstitch/bidirs/hist/", mode: 'copy', pattern: "*.html"
    publishDir "${params.workdir}/fstitch/bidirs/stats/", mode: 'copy', pattern: "*.txt"

    when:
    params.fstitch

    input:
    set val(prefix), file(bg) from fstitch_bg
    set val(prefix), file(pos_bg) from pos_fstitch_bg
    set val(prefix), file(neg_bg) from neg_fstitch_bg

    output:
    file ("*.hmminfo") into fs_train_out
    file ("*.fstitch_seg.bed") into fs_seg_out
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
        -o ${prefix}.pos.fstitch_seg.bed

    ${params.fstitch_path} segment \
        -s - \
        -b ${neg_bg} \
        -p ${prefix}.fstitch.hmminfo \
        -o ${prefix}.neg.fstitch_seg.bed

    cat ${prefix}.pos.fstitch_seg.bed \
        ${prefix}.neg.fstitch_seg.bed \
        | sortBed > ${prefix}.fstitch_seg.bed

    bidir \
        -b ${prefix}.fstitch_seg.bed \
        -g ${params.genome_refseq} \
        -o ${prefix}.fstitch_bidir.bed \
        -p \
        -s
    """
}

println "[Log 3]: Done Running FStitch\n"


// PART 4: Running TFit

process tfit {
    println "[Log 4]: Running TFit full model"

    tag "$prefix"
    memory '25 GB'
    time '20h'
    cpus 16
    queue 'short'
    validExitStatus 0

    publishDir "${params.workdir}/tfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.workdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"

    when:
    params.tfit

    input:
    set val(prefix), file (bidirs) from fs_bidir_out
    set val(prefix), file(bg) from tfit_bg

    output:
    set val(prefix), file ("${prefix}.tfit_bidirs.bed") into tfit_bed_out
    file ("*.tsv") into tfit_full_model_out
    file ("*.log") into tfit_logs_out

    script:
        """
        export OMP_NUM_THREADS=16

        ${params.tfit_path} model \
            -bg ${bg} \
            -s ${bidirs} \
            -N $prefix \
            -l \
            -o ${prefix}.tfit_bidirs.bed \
            --threads 16 \
        """
}



process prelimtfit {
    println "[Log 4]: Running TFit prelim"

    tag "$prefix"
    memory '100 GB'
    time '48h'
    cpus 16
    queue 'long'
    validExitStatus 0

    publishDir "${params.workdir}/prelimtfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.workdir}/prelimtfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    publishDir "${params.workdir}/prelimtfit/prelim", mode: 'copy', pattern: "*tfit_prelim.bed"

    when:
    params.prelimtfit

    input:
    set val(prefix), file(bg) from prelimtfit_bg

    output:
    file ("*tfit_prelim.bed") into tfit_prelim_out
    file ("*tfit_bidirs.bed") into prelimtfit_bed_out
    file ("*.tsv") into prelimtfit_full_model_out
    file ("*.log") into prelimtfit_logs_out

    script:
        """
        export OMP_NUM_THREADS=16

        ${params.tfit_path} prelim \
            -bg ${bg} \
            -N $prefix \
            -l \
            -o ${prefix}.tfit_prelim.bed \
            --threads 16

        ${params.tfit_path} model \
            -bg ${bg} \
            -s ${prefix}.tfit_prelim.bed \
            -N $prefix \
            -l \
            -o ${prefix}.tfit_bidirs.bed \
            --threads 16
        """
}

println "[Log 4]: Done Running TFit\n"



// PART 5: Preparing bigwig files for dREG

process dreg_prep {
    println "[Log 5]: Generating bigwig files for dREG"

    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$prefix"
    memory '150 GB'
    queue 'short'

    publishDir "${params.workdir}/bigwig/", mode: 'copy', pattern: "*.bw"

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
    memory '120 GB'
    time '14h'
    cpus 4
    queue 'titan'
    validExitStatus 0

    publishDir "${params.workdir}/dreg/", mode: 'copy', pattern: "*dREG*"

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
