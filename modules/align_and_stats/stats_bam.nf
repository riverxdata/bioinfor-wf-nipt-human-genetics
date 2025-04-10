process STATS_BAM {
    label 'process_low'
    tag "$sample_id"
    container 'quay.io/biocontainers/samtools:1.15.1--h6899075_1'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.bamstats")

    script:
    """
    samtools stats ${bam} > ${sample_id}.sorted.rmdup.realign.BQSR.bamstats
    """
}
