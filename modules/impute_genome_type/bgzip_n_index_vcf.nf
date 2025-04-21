process BGZIP_AND_INDEX_VCF {
    label 'process_low'
    tag "$chr"
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(chr), path(vcf)

    output:

    tuple val(chr), path("${vcf}.gz"), path("${vcf}.gz.csi")

    script:
    """
    bgzip ${vcf}
    bcftools index -f ${vcf}.gz
    """
}
