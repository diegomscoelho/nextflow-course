process KALLISTO_INDEX {

    input:
    path transcriptome

    output:
    path 'kallisto_index'

    script:
    """
    kallisto index -i kallisto_index $transcriptome
    """
}

process KALLISTO_QUANT {
    tag "Kallisto on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path kallisto_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    kallisto quant -t $task.cpus -i $kallisto_index  -o $sample_id ${reads[0]} ${reads[1]}
    """
}
