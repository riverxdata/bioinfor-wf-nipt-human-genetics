process BEDTOOLS_GENOMECOV {
    label 'process_medium'
    tag "$sample_id"
    container 'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz"), path("${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz.tbi")

    script:
    """
    bedtools genomecov -ibam ${bam} -bga -split | bgzip > ${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz
    tabix -p bed ${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz
    """
}
