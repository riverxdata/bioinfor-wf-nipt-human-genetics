process INDEX_BED {
    label 'process_medium'
    tag "$sample_id"
    container 'quay.io/biocontainers/tabix:0.2.5--0'

    input:
    tuple val(sample_id), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.{gz,gz.tbi}")

    script:
    """
    tabix -p bed $bed
    """
}
