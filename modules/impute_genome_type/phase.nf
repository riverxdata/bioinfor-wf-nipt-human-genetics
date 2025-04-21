process PHASE {
    label 'process_low'
    container 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521'
    tag "${merged_vcf[0].getName()}"

    input:
    tuple val(chr), path(gmap), path(merged_vcf), path(ref_vcf), path(chunk_txt)

    output:
    tuple val(chr), path("imputed.glimpse.${chunk_txt.baseName}.vcf")

    script:
    """
    IRG=\$(cut -f3 ${chunk_txt})
    ORG=\$(cut -f4 ${chunk_txt})
    /usr/bin/GLIMPSE_phase_v1.1.1 \\
        --input ${merged_vcf[0]} \\
        --reference ${ref_vcf[0]} \\
        --map ${gmap} \\
        --input-region \$IRG \\
        --output-region \$ORG \\
        --output imputed.glimpse.${chunk_txt.baseName}.vcf
    """
}
