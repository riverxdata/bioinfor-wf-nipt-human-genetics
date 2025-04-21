process PREPARE_REF_VCF{
    label 'process_low'
    tag "$vcf"
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
    tuple val(chr), path(vcf), path(vcf_index)

    output:
    tuple val(chr), path("biallelic.snp.maf0.001.${chr}.vcf.gz*"), emit: vcf
    tuple val(chr), path("biallelic.snp.maf0.001.${chr}.sites.tsv.gz"), emit: tsv
    tuple val(chr), path("biallelic.snp.maf0.001.${chr}.vcf.gz*"), path("biallelic.snp.maf0.001.${chr}.sites.tsv.gz"), emit: all

    script:
    """
    # Conduct normalization and filtration of the reference panel
    bcftools norm -m -any ${vcf} -Ou --threads ${task.cpus} | \
    bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads ${task.cpus} -Oz -o biallelic.snp.maf0.001.${chr}.vcf.gz
    bcftools index -f biallelic.snp.maf0.001.${chr}.vcf.gz


    # Extracting variable positions in the reference panel
    bcftools view -G -m 2 -M 2 -v snps biallelic.snp.maf0.001.${chr}.vcf.gz -Oz -o biallelic.snp.maf0.001.${chr}.sites.vcf.gz --threads ${task.cpus}
    bcftools index -f biallelic.snp.maf0.001.${chr}.sites.vcf.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' biallelic.snp.maf0.001.${chr}.sites.vcf.gz | bgzip -c > biallelic.snp.maf0.001.${chr}.sites.tsv.gz
    tabix -s1 -b2 -e2 biallelic.snp.maf0.001.${chr}.sites.tsv.gz
    """
}