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



// Configurable variables
params.name = false
params.email = false
params.genome = 
params.plaintext_email = false
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
//params.workdir = "./nextflowTemp"
//output_docs = file("$baseDir/docs/output.md")

////////////////////////////////////////////////////
/* --          INITIALIZE LOG                  -- */
////////////////////////////////////////////////////
log.info "====================================="
log.info "  Bidirectional Detection Pipeline   "
log.info "====================================="
log.info "name         : ${params.name}"
log.info "email        : ${params.email}"
log.info "\n"


// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }


if ( params.fstitch_path ){
    fstitch_path = file("${params.fstitch_path}")
}

if ( params.fstitch_train){
    fstitch_train = file("${params.fstitch_train}")
}

if ( params.tfit_path){
    tfit_path = file("${params.tfit_path}")
}

if ( params.dreg_path){
   dreg_path = filr("${params.dreg_path}")
}


/*
 *STEP 1 - Create bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    memory '40 GB'
    time '4h'
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"      

    input:
    set val(name), file(bam_file) from sorted_bam_for_bedgraph
    set val(name), file(bam_indices) from sorted_bam_index_for_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("${name}.pos.bedGraph") into pos_non_normalized_bedgraphs, pos_fstitch_bg
    set val(name), file("${name}.neg.bedGraph") into neg_non_normalized_bedgraphs, neg_fstitch_bg
    set val(name), file("${name}.bedGraph") into non_normalized_bedgraphs, fstitch_bg, tfit_bg, prelimtfit_bg, nqc_bg
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    if (params.singleEnd) {
    """    
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        > ${name}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.neg.bedGraph
    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph
    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph
    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    """
    } else {
    """   
    samtools view \
        -h -b -f 0x0040 \
        ${bam_file} \
        > ${name}.first_pair.bam
        
    samtools view \
        -h -b -f 0x0080 \
        ${bam_file} \
        > ${name}.second_pair.bam
        
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        > ${name}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.first_pair.neg.bedGraph
                     
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | sortBed \
        > ${name}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${name}.second_pair.neg.bedGraph
                     
    unionBedGraphs \
        -i ${name}.first_pair.pos.bedGraph ${name}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.pos.bedGraph 
        
    unionBedGraphs \
        -i ${name}.first_pair.neg.bedGraph ${name}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.neg.bedGraph 
    
    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph
    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph
    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph     
    """
    }
 }

/*
 *STEP 2 - Create bedGraphs and bigwigs for dREG
 */

process dreg_prep {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    memory '150 GB'
    publishDir "${params.outdir}/mapped/dreg_input", mode: 'copy', pattern: "*.bw"
    
    when:
    params.dreg

    input:
    set val(name), file(cram_file) from cram_dreg_prep
    set val(name), file(cram_index) from cram_index_dreg_prep
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.bw") into dreg_bigwig

    script:
    """
    echo "Creating BigWigs suitable as inputs to dREG"
    
    export CRAM_REFERENCE=${genome}    
    
    bamToBed -i ${cram_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sortBed -i ${name}.dreg.bed > ${name}.dreg.sort.bed
    
    echo positive strand processed to bedGraph
    
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand + \
            > ${name}.pos.bedGraph
    sortBed \
            -i ${name}.pos.bedGraph \
            > ${name}.pos.sort.bedGraph
            
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${name}.neg.bedGraph
    sortBed \
            -i ${name}.neg.bedGraph \
            > ${name}.neg.sort.bedGraph
            
    echo negative strand processed to bedGraph
    
    ${params.bedGraphToBigWig} ${name}.pos.sort.bedGraph ${chrom_sizes} ${name}.pos.bw
    ${params.bedGraphToBigWig} ${name}.neg.sort.bedGraph ${chrom_sizes} ${name}.neg.bw
    
    echo bedGraph to bigwig done
    """
 }



/*
 * STEP 3 - Running FStitch
 */

process FStitch {
    tag "$name"
    memory '50 GB'
    time '1h'
    validExitStatus 0
    publishDir "${params.outdir}/fstitch/", mode: 'copy', pattern: "*.hmminfo"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "*.fstitch_seg.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "${name}.fstitch_bidir.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "*fstitch_bidir.{short,long}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/hist/", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/fstitch/bidirs/stats/", mode: 'copy', pattern: "*.txt"
    
    when:
    params.fstitch
    
    input:
    set val(name), file(bg) from fstitch_bg
    set val(name), file(pos_bg) from pos_fstitch_bg
    set val(name), file(neg_bg) from neg_fstitch_bg
        
    output:
    file ("*.hmminfo") into fs_train_out
    file ("*.fstitch_seg.bed") into fs_seg_out
    set val(name), file ("*.fstitch_bidir.bed") into fs_bidir_out
    file ("*fstitch_bidir.{short,long}.bed") into fs_bidir_short_long_out
    file ("*.html") into fs_bidir_plot_out
    file ("*.txt") into fs_bidir_stats_out
    
    script:
    """
    ${fstitch_path} train \
        -s + \
        -b ${bg} \
        -t ${fstitch_train} \
        -o ${name}.fstitch.hmminfo
        
     ${fstitch_path} segment \
        -s + \
        -b ${pos_bg} \
        -p ${name}.fstitch.hmminfo \
        -o ${name}.pos.fstitch_seg.bed
     ${fstitch_path} segment \
        -s - \
        -b ${neg_bg} \
        -p ${name}.fstitch.hmminfo \
        -o ${name}.neg.fstitch_seg.bed
    cat ${name}.pos.fstitch_seg.bed \
        ${name}.neg.fstitch_seg.bed \
        | sortBed > ${name}.fstitch_seg.bed
    bidir \
        -b ${name}.fstitch_seg.bed \
        -g ${genome_refseq} \
        -o ${name}.fstitch_bidir.bed \
        -p \
        -s    
    """
}


/*
 * STEP 4 - Running Tfit
 */

process tfit {
    tag "$name"
    memory '25 GB'
    time '20h'
    cpus 16
    queue 'short'
    validExitStatus 0
    publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    
    when:
    params.tfit
    
    input:
    set val(name), file (bidirs) from fs_bidir_out
    set val(name), file(bg) from tfit_bg
        
    output:
    set val(name), file ("${name}.tfit_bidirs.bed") into tfit_bed_out
    file ("*.tsv") into tfit_full_model_out
    file ("*.log") into tfit_logs_out
        
    script:
        """
        export OMP_NUM_THREADS=16
        
        ${tfit_path} model \
            -bg ${bg} \
            -s ${bidirs} \
            -N $name \
            -l \
            -o ${name}.tfit_bidirs.bed \
            --threads 16 \
        """
}

process prelimtfit {
    tag "$name"
    memory '100 GB'
    time '48h'
    cpus 16
    queue 'long'
    validExitStatus 0
    publishDir "${params.outdir}/prelimtfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.outdir}/prelimtfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    publishDir "${params.outdir}/prelimtfit/prelim", mode: 'copy', pattern: "*tfit_prelim.bed"
    
    when:
    params.prelimtfit
    
    input:
    set val(name), file(bg) from prelimtfit_bg
        
    output:
    file ("*tfit_prelim.bed") into tfit_prelim_out
    file ("*tfit_bidirs.bed") into prelimtfit_bed_out
    file ("*.tsv") into prelimtfit_full_model_out
    file ("*.log") into prelimtfit_logs_out
        
    script:
        """
        export OMP_NUM_THREADS=16
        
        ${tfit_path} prelim \
            -bg ${bg} \
            -N $name \
            -l \
            -o ${name}.tfit_prelim.bed \
            --threads 16
        
        ${tfit_path} model \
            -bg ${bg} \
            -s ${name}.tfit_prelim.bed \
            -N $name \
            -l \
            -o ${name}.tfit_bidirs.bed \
            --threads 16
        """
}