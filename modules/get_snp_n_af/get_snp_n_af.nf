process GET_SNP_AND_AF {
    container 'docker.io/nttg8100/basevar:1.2.2'
    cpus 32
    memory "200.GB"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path(fasta_index)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")
    tuple val(sample_id), path("${sample_id}.cvg.tsv.gz")

    script:
    """
    REF=\$(find -L ./ -name "*.fasta") 
    /opt/basevar/bin/basevar basetype -R \$REF \
        -B 200 \
        -t ${task.cpus} \
        -I ${bam[0]} \
        -r chr22 \
        --output-vcf ${sample_id}.vcf.gz \
        --output-cvg ${sample_id}.cvg.tsv.gz
    """
}
