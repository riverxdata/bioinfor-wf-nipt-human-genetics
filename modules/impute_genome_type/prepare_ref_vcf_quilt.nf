process PREPARE_REF_VCF_QUILT{
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

    zcat biallelic.snp.maf0.001.${chr}.vcf.gz| grep '^#'> biallelic.snp.maf0.001.unique_pos.vcf 
    zcat biallelic.snp.maf0.001.${chr}.vcf.gz| awk '!/^#/ {print \$0}'|awk 'length(\$4)==1'|awk 'length(\$5)==1'| awk -F '\\t' '!a[\$2]++' >> biallelic.snp.maf0.001.unique_pos.vcf  
    bgzip -f biallelic.snp.maf0.001.unique_pos.vcf

    # convert format
    bcftools convert --haplegendsample biallelic.snp.maf0.001.unique_pos biallelic.snp.maf0.001.unique_pos.vcf.vcf.gz 
    sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' biallelic.snp.maf0.001.unique_pos.samples
    """
}