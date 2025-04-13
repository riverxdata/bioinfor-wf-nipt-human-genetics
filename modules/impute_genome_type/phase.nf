$GLIMPSE_phase --input ${VCF} --reference ${REF} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
$bgzip ${OUT}
$bcftools index -f ${OUT}.gz
process PHASE {
    label 'process_medium'
    container 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521'
    tag "${vcf[0]}"

    input:
    tuple val(chr), path(vcf)

    output:
    tuple val(chr), path("chunks.G10K.${chr}.txt")

    script:
    """
    /usr/bin/GLIMPSE_phase_v1.1.1 \\
        --input ${vcf[0]} \\
        --reference \${REF} \\
        --map ${MAP} \\
        
    """
}