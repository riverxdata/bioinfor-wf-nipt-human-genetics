process GET_ACCURACY {
    label 'process_low'
    tag "$chr"
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(chr), path(vcf)

    output:

    tuple val(chr), path("${vcf}.gz"), path("${vcf}.gz.csi")

    script:
    """
    bcftools stats \
        $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.txt
    """
}
