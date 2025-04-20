process LIGATE {
    label 'process_medium'
    container 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521'
    tag "${chr}"

    input: 
    tuple val(chr), path(vcf_list), path(vcf_index_list)

    output:
    tuple val(chr), path("high_dep_100.${chr}_imputed.vcf")
    script:
    """
    ls imputed.glimpse.*.vcf.gz > high_dep_100.${chr}_imputed_list.txt
    /usr/bin/GLIMPSE_ligate_v1.1.1 --input high_dep_100.${chr}_imputed_list.txt --output high_dep_100.${chr}_imputed.vcf
    """
}
