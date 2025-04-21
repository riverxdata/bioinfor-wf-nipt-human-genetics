process PLINK2_PCA {
    label 'process_medium'
    tag "all"
    container 'quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_4'

    input:
    path(vcf_and_index)

    output:
    path("plink2_genotype*")

    script:
    """
    # Step 1: Convert VCF to PLINK binary format
    plink2 \
        --vcf *.vcf.gz \
        --make-bed \
        --out plink2_genotype \
        --threads ${task.cpus}

    # Step 2: LD-pruning
    plink2 \
        --bfile plink2_genotype \
        --maf 0.01 \
        --threads ${task.cpus} \
        --indep-pairwise 500 50 0.2 \
        --out plink2_ld_pruned \
        --bad-ld

    # Step 3: Remove related samples using KING cutoff
    plink2 \
        --bfile plink2_genotype \
        --extract plink2_ld_pruned.prune.in \
        --king-cutoff 0.177 \
        --out plink2_remove_related_samples \
        --threads ${task.cpus}

    # Step 4: Get reference allele frequencies
    plink2 \
        --bfile plink2_genotype \
        --freq \
        --out plink2_genotype \
        --threads ${task.cpus}

    # Step 5: PCA using unrelated and LD-pruned samples
    plink2 \
        --bfile plink2_genotype \
        --keep plink2_remove_related_samples.king.cutoff.in.id \
        --extract plink2_ld_pruned.prune.in \
        --read-freq plink2_genotype.afreq \
        --threads ${task.cpus} \
        --pca approx allele-wts 3 \
        --out plink2_genotype

    # Step 6: Project all samples onto PCs
    plink2 \
        --bfile plink2_genotype \
        --threads ${task.cpus} \
        --read-freq plink2_genotype.afreq \
        --score plink2_genotype.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
        --score-col-nums 6-8 \
        --out plink2_genotype
    """
}
