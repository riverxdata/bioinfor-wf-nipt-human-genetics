process INDEX_BAM {
    label 'process_low'
    tag "$sample_id"
    container 'quay.io/biocontainers/samtools:1.15.1--h6899075_1'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}
