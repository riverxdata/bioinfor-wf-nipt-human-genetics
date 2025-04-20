process PLINK2_GWAS{
    label 'process_low'
    tag "all"
    container 'quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_4'

    input:
    path(phenotype)
    path(pca)

    output:
    path("plink2_gwas*")

    script:
    """
    plink2 \
        --bfile ${pca[0].baseName} \
        --pheno ${phenotype} \
        --pheno-name B1 \
        --maf 0.01 \
        --covar ${pca[0].baseName}.sscore \
        --covar-col-nums 5-7 \
        --glm hide-covar firth firth-residualize single-prec-cc \
        --threads ${task.cpus} \
        --out plink2_gwas
    """
}
