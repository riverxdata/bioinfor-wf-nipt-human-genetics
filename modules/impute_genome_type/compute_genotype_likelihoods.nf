process COMPUTE_GENOTYPE_LIKELIHOODS{
    label 'process_medium'
    tag "${sample_id}-${chr}"
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    tuple val(chr), path(vcf)
    tuple val(chr), path(tsv)
    path ref

    output:
    tuple val("${sample_id}"), val("${chr}"),path("genotype_likelihoods_${sample_id}.${chr}.vcf.gz"), path("genotype_likelihoods_${sample_id}.${chr}.vcf.gz.csi")

    script:
    """
    REF=\$(find -L ./ -name "*.fasta") 
    bcftools mpileup -f \${REF} -I -E -a 'FORMAT/DP' -T ${vcf[0]} -r ${chr} ${bam} -Ou | bcftools call -Aim -C alleles -T ${tsv} -Oz -o  genotype_likelihoods_${sample_id}.${chr}.vcf.gz
    bcftools index -f genotype_likelihoods_${sample_id}.${chr}.vcf.gz
    """
}