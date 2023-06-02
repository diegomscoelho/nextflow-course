/*
 * pipeline input parameters
 */
params.reads = "$projectDir/samples/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/samples/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */

include { INDEX; QUANTIFICATION } from "$projectDir/nf_modules/salmon.nf"
include { KALLISTO_INDEX; KALLISTO_QUANT } from "$projectDir/nf_modules/kallisto.nf"
include { FASTQC; MULTIQC } from "$projectDir/nf_modules/qc.nf"

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    kallisto_index_ch = KALLISTO_INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    kallisto_quant_ch = KALLISTO_QUANT(kallisto_index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
