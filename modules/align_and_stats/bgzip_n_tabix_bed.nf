process BGZIP_AND_TABIX_BED {
    label 'process_low'
    tag "$sample_id"
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(sample_id), path(bed)

    output:

    tuple val(sample_id), path("${bed}.gz"), path("${bed}.gz.tbi")

    script:
    """
    bgzip ${bed}
    tabix -p bed ${bed}.gz
    """
}
