#!/usr/bin/env nextflow
/*
========================================================================================
                         Long read mapping
========================================================================================
 Long read mapping pipeline. This read alignment requires a fastq read file and a fasta reference genome 

 Colin Davenport 

 #### Homepage / Documentation Changelog

v0.1.1  
v0.1.0  Generate stats and multiqc
v0.0.9  Change output names, allow only one fastq input (else waits during collect steps after longshot)
v0.0.8  bcf sorting
v0.0.7  Adding bctftools concat for gathering VCF, and bam to cram conversion
v0.0.6  Longshot scatter to vcf with flatten and each working
v0.0.5  Trial -r region variant calling acceleration (alternative: bedtools makewindows on .fai)
v0.0.4  Add memory scaling cluster retries 
v0.0.3  Add longshot SNV calling
v0.0.2  Minimap2 and samtools working
v0.0.1  First working version


// TODO - following short read pipeline
// samtools stats --threads {threads} {input.bam} > {output.stats}
// multiqc         multiqc --force --no-data-dir -n {output.html} .
// vt normalize {input.bcf} -r {config[reference]} | vt uniq - -o {output.bcf} 2> {log}   // but is this sensible for long read SNVs?
// filters hard = not pass:         bcftools view --apply-filters PASS -Ob --threads {threads} -o {output.bcf} {input.bcf}
// soft = v complex, strand bias, etc

----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is:
      conda activate nextflow
      nextflow run longread_alignment.nf  --fasta /path/to/x.fa 
      //nextflow run longread_alignment --genome ZR_S706_v1 -profile conda 

    Arguments:
      --fasta [file]                  Path to Fasta reference. (Default: false)
      --fai [file]                    Path to Fasta.fai reference index used for regions specific SNV calling. (Default: false)
      --fastq.gz [str]                fastq.gz read set
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, singularity

    Other
      --outdir [file]                 The output directory where the results will be saved (Default: './output')
    """.stripIndent()
}

// use modern nextflow
nextflow.enable.dsl = 2

/*
* Parameter defaults
*/

params.help=false
params.save_align_intermeds=true
params.outdir = "output"
params.publish_dir_mode = "copy"
params.fastq = "/mnt/beegfs/scratch/bioinformatics/colin/dev/longread_alignment_nf/test_beet_big.fastq.gz"
params.fasta = "/mnt/omics/data/projects/26890/uploaded/214815-ZR_S706_v1_chromosomes.fasta"
params.fai = "/mnt/omics/data/projects/26890/uploaded/214815-ZR_S706_v1_chromosomes.fasta.fai"
params.test = "yes"
params.min_cov = ""
params.min_alt_count = ""
params.min_alt_frac = ""


if (params.test.toString().contains("yes")) {
    println("\nSNP calling defaults set very low for testing only!\n") 
    params.min_cov = 1
    params.min_alt_count = 1
    params.min_alt_frac = 0.25
}
else if (params.test.contains("no")) {
    println("\nSNP calling defaults set to normal for production!\n") 
    params.min_cov = 10
    params.min_alt_count = 6
    params.min_alt_frac = 0.25
}
else {
    print("\nWarning: params.test setting failed !\n")
}

//script:
//    params.fai = params.fasta + ".fai"

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
* Workflow
*/


workflow {

    println "Starting longread_alignment.nf"
    println "Version 0.1.0 by Colin Davenport"

    // File inputs
    input_fastq = Channel.fromPath(params.fastq, checkIfExists: true)

    chunksize = Channel.value(1000)


    // run processes
    
    //prepare_fai(params.fai)

    // Run SNP calling per chromosome (might be easier to split bams, then per bam ?)
    // Alternative - use per window: bedtools makewindows -g 214815-ZR_S706_v1_chromosomes.fasta.fai -w 10000000 > windows.fai
    // limits to .take(x) chromosomes for dev only
    println "\n ###### Warning ! Will only perform SNP calling for these chromosomes: ######"
    fai_channel = Channel
        .fromPath(params.fai)
        .splitText( each:{ it.split().first() } )
        .take(9)
        .map { [ it ] }
        .view()
    
    minimap2_align(input_fastq)
    sort_bam(minimap2_align.out)

    // generate alignment stats
    bam_stats(sort_bam.out)

    // convert bam to cram format
    convert_bam_cram(sort_bam.out)
    
    // longshot, run split by chromosome -from bam
    longshot(sort_bam.out, fai_channel.flatten())
    
    // longshot, run split by chromosome - from cram -smaller, preferred
    //longshot(convert_bam_cram.out, fai_channel.flatten())
 
    // concatenate vcf.gz to bcf. 
    //concat_bcf(reheader_vcfs.out)
    concat_bcf(longshot.out.collect())

    // sort bcf results
    sort_bcf(concat_bcf.out)

    // generate bcf stats
    stats(sort_bcf.out)

    // multiqc
    multiqc(stats.out.collect(), bam_stats.out)


}




/*
 *  Run minimap2
 */

process minimap2_align {

    cpus = 16
	// If job fails, try again with more memory
	memory { 32.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/wochenende/'


    tag "$name"
    label 'process_medium'
    
    // Do not save sam files to output folder by default
    /*
    *if (params.save_align_intermeds) {
    *    publishDir path: "${params.outdir}/minimap2", mode: params.publish_dir_mode,
    *        saveAs: { filename ->
    *                      if (filename.endsWith('.flagstat')) "$filename"
    *                      else filename
    *                }
    *}
    */


    input:
    file fastq


    output:
    file "${prefix}.sam"
    

    script:
    prefix = fastq.name.toString().tokenize('.').get(0)
    name = fastq

    """
    minimap2 -x map-ont -a --split-prefix ${prefix} -t $task.cpus  -o ${prefix}.sam $params.fasta $fastq
    """
}


/*
 *  Convert BAM to coordinate sorted BAM, make stats, flagstat, idxstats
 */

process sort_bam {

    cpus = 8
	// If job fails, try again with more memory
	memory { 32.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.flagstat')) "$filename" 
                          else if (filename.endsWith('.idxstats')) "$filename" 
                          else if (filename.endsWith('.stats')) "$filename" 
                          else filename
                    }
    }


    input:
    file bam


    output:
    file "${prefix}.sorted.bam"
    file "${prefix}.sorted.bam.bai"
    file "${prefix}.sorted.bam.flagstat"
    file "${prefix}.sorted.bam.idxstats"
    file "${prefix}.sorted.bam.stats"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}



/*
 *  make bam stats, flagstat, idxstats
 */

process bam_stats {

    cpus = 8
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.flagstat')) "$filename" 
                          else if (filename.endsWith('.idxstats')) "$filename" 
                          else if (filename.endsWith('.stats')) "$filename" 
                          else filename
                    }
    }


    input:
    file bam


    output:
    file "${prefix}.sorted.bam.flagstat"
    file "${prefix}.sorted.bam.idxstats"
    file "${prefix}.sorted.bam.stats"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}


/*
 *  Convert BAM to coordinate sorted BAM, make stats, flagstat, idxstats
 */

process convert_bam_cram {

    cpus = 8
	// If job fails, try again with more memory
	memory { 32.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/wochenende/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.cram')) "$filename"                         
                          else filename
                    }
    }


    input:
    file bam
    file bai
    file flagstat
    file idxstats
    file stats

    output:
    file "${prefix}.cram"
    file "${prefix}.cram.crai"


    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    // cram conversion  samtools view -@ {threads} -C -T {params.reference} -o {output.cram} {output.bam}

    """
    samtools view -@ $task.cpus -C -T $params.fasta -o ${prefix}.cram $bam
    samtools index ${prefix}.cram
    """
}




/*
 *  Call SNVs from long read BAM
 *  TODO - parallelize by chr, then recombine?
 */

process longshot {
    //echo true
    cpus = 2
	// If job fails, try again with more memory
	memory { 16.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name + $region"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/longshot", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.vcf.gz')) "$filename"
                          else if (filename.endsWith('.tbi')) "$filename"
                          else filename
                    }
    }


    input:
    file sorted_bam
    file bai
    each region
  
    output:
    //file vcf_gz
    file vcf_reheader_gz

    script:
    prefix = sorted_bam.toString().tokenize('.').get(0)
    name = sorted_bam
    vcf =  prefix + "_" + region + ".vcf"
    vcf_gz = vcf + ".gz"
    vcf_reheader_gz = vcf + "_reheader.gz"
    //print(params.fai)
    
    """
    longshot --region ${ region } --bam $sorted_bam --ref $params.fasta  --min_cov $params.min_cov --min_alt_count $params.min_alt_count --min_alt_frac $params.min_alt_frac --sample_id ${prefix} --out ${vcf} > ${prefix}${region + "_out.txt"} 2> ${prefix}${region + "_err.txt"}
    
    bgzip -f -k -@ $task.cpus ${vcf}
    tabix ${vcf + ".gz"}
    
    bcftools reheader -f ${params.fai} -o $vcf_reheader_gz   ${vcf + ".gz"}

    
    """    

    //longshot --region ${ region } --bam $sorted_bam --ref $params.fasta --out ${prefix + "_"}${region}.vcf --min_cov $params.min_cov --min_alt_count $params.min_alt_count --min_alt_frac $params.min_alt_frac --sample_id ${prefix} > ${prefix}${region + "_out.txt"} 2> ${prefix}${region + "_err.txt"}

    // forward, reverse, merge agreement as D. Eccles and Janinas suggestion?
    // or set strand_bias_pvalue_cutoff higher ?
    // longshot --anchor_length <int> --band_width <Band width> --density_params <string> --hap_converge_delta <float> --hap_assignment_qual <float> --het_snv_rate <float> --hom_snv_rate <float> --bam <BAM> --ref <FASTA> --max_cigar_indel <int> --max_cov <int> --max_window <int> --min_allele_qual <float> --min_cov <int> --min_mapq <int> --out <VCF> --potential_snv_cutoff <float> --min_alt_count <int> --min_alt_frac <float> --sample_id <string> --strand_bias_pvalue_cutoff <float> --max_snvs <int> --ts_tv_ratio <float>
    // longshot defaults are permissive, suggest more restrictive as with short read pipeline! eg 10, 6, 0.3 respectively ?
    // --min_cov = 6
    // -e, --min_alt_count <int>                  Require a potential SNV to have at least this many alternate allele observations. [default: 3]
    // -E, --min_alt_frac <float>                 Require a potential SNV to have at least this fraction of alternate allele observations. [default: 0.125]
    


}



/*
 *  Combine VCFs, reheader, create bcf from long read BAM
 */

process concat_bcf {

    cpus = 8
	// If job fails, try again with more memory
	memory { 16.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/bcftools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.bcfxxxx')) "$filename"
                          else filename
                    }
    }

    input:
    file vcflist
    
    
    output:
    file bcf
    file bcfindex

    script:
    //prefix = vcflist[0].toString().tokenize('.').get(0)
    name = "concatenate"
    first_vcf = vcflist[0]
    bcf = "concatenated.bcf"
    bcfindex = "concatenated.bcf.csi"
    sort_dir = "sort_tmp"

    """
    bcftools concat --threads $task.cpus -Ov -o concatenated.vcf.gz  ${vcflist}
    tabix concatenated.vcf.gz

    bcftools reheader -f ${params.fai} -o concatenated2.vcf.gz concatenated.vcf.gz
    tabix concatenated2.vcf.gz

    bcftools view -Ob -o concatenated.bcf concatenated2.vcf.gz
    bcftools index concatenated.bcf
    
    """


}




/*
 *  Sort bcf
 */

process sort_bcf {

    cpus = 2
	// If job fails, try again with more memory
	memory { 16.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/bcftools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.s.bcf')) "$filename"
                          else if (filename.endsWith('.s.csi')) "$filename"
                          else filename
                    }
    }

    input:
    file bcf
    file bcf_index
    
    
    output:
    file bcf_s
    file bcf_s_index

    script:
    prefix = params.fastq.toString().tokenize('.').get(0)
    name = prefix
    bcf_s = prefix + ".s.bcf"
    bcf_s_index = prefix + ".s.bcf.csi"
    sort_dir = "sort_tmp"

    """

    mkdir -p ${sort_dir}
    bcftools sort --temp-dir ${sort_dir} -Ob -o ${prefix + ".s.bcf"} concatenated.bcf
    bcftools index ${prefix + ".s.bcf"}
    """



}



/*
 *  Generate stats
 */

process stats {
    cpus = 8
	// If job fails, try again with more memory
	memory { 4.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/bcftools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.txt')) "$filename"
                          else filename
                    }
    }
    input:
    file bcf
    file bcf_index
    
    output:
    path("*.txt")


    script:
    prefix = bcf.toString().tokenize('.').get(0)
    name = prefix


    """
    vt peek $bcf > ${prefix + "vtpeek.txt"}
    bcftools stats --threads $task.cpus $bcf > ${prefix + "bcftools_stats.txt"}

    """
    
}


/*
 *  Multiqc
 */

process multiqc {
    cpus = 2
	// If job fails, try again with more memory
	memory { 4.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/multiqc", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.html')) "$filename"
                          else filename
                    }
    }

    input:
    path multiqc_files
    file flagstat
    file idxstats
    file stats

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    //when:
    //task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    name = "All stats files"
    
    """
    multiqc -f .
    
    """

}
