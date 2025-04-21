process CONCAT_VCF {
    tag "$prefix"
    label 'process_high'
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    val(prefix)
    path vcf

    output:
    path("${prefix}.vcf.gz*")

    script:
    """
    bcftools concat --threads ${task.cpus} -a -Oz -o ${prefix}.vcf.gz *.vcf.gz
    bcftools index --threads ${task.cpus} ${prefix}.vcf.gz
    """
}
