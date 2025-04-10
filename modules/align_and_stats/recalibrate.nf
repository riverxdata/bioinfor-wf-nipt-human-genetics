process RECABLIRATOR {
    label 'process_medium'
    tag "$sample_id"
    container 'quay.io/biocontainers/gatk:3.8--0'

    input:
    tuple val(sample_id), path(bam), path(reference), path(dbsnp), path(known_indels1), path(known_indels2)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.bam")

    script:
    """
    gatk_jar="/usr/local/bin/gatk.jar"

    java -Xmx${task.memory.toMega()}m -jar \$gatk_jar \
        -T BaseRecalibrator \
        -nct ${task.cpus} \
        -R ${reference} \
        -I ${bam} \
        --knownSites ${dbsnp} \
        --knownSites ${known_indels1} \
        --knownSites ${known_indels2} \
        -o ${sample_id}.recal_data.table

    java -Xmx${task.memory.toMega()}m -jar \$gatk_jar \
        -T PrintReads \
        -nct ${task.cpus} \
        -R ${reference} \
        --BQSR ${sample_id}.recal_data.table \
        -I ${bam} \
        -o ${sample_id}.sorted.rmdup.realign.BQSR.bam
    """
}
