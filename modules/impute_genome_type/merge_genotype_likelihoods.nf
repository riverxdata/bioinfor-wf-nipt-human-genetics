process MERGE_GENOTYPE_LIKELIHOODS{
    tag "$chr"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(chr), path(vcf), path(vcf_index)

    output:
    tuple val(chr),path("high_dep_100.${chr}.vcf.gz*")

    script:
    """
    ls *.vcf.gz > high_dep_100.chr${chr}_GL_list.txt
    bcftools merge --force-single -m none -r ${chr} -Oz -o high_dep_100.${chr}.vcf.gz -l high_dep_100.chr${chr}_GL_list.txt
    bcftools index -f high_dep_100.${chr}.vcf.gz
    """
}