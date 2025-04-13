process BEDTOOLS_GENOMECOV {
    label 'process_medium'
    tag "$sample_id"
    container 'docker.io/pegi3s/bedtools'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed")

    script:
    """
    bedtools genomecov -ibam ${bam} -bga -split > ${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed
    """
}
