process GET_IMPUTATION_CHUNK {
    label 'process_medium'
    container 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521'
    tag "${vcf[0]}"

    input:
    tuple val(chr), path(vcf)

    output:
    tuple val(chr), path("region_chunk_${chr}_*.txt")

    script:
    """
    /usr/bin/GLIMPSE_chunk_v1.1.1 \\
        --input ${vcf[0]} \\
        --region ${chr} \\
        --window-size 2000000 \\
        --buffer-size 200000 \\
        --output chunks.G10K.${chr}.txt
    awk '{print > ("region_chunk_${chr}_" \$1 ".txt")}' chunks.G10K.${chr}.txt
    """
}