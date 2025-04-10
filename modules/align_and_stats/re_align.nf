process REALIGNER {
    label 'process_medium'
    tag "$sample_id"
    container 'quay.io/biocontainers/gatk:3.8_0--hdfd78af_11'

    input:
    tuple val(sample_id), path(bam)
    path reference
    path known_indels_1
    path known_indels_2

    output:
    path "${sample_id}.indel_target_intervals.list"

    script:
    """
    java -Xmx${task.memory.toMega()}m -jar \$GATK_HOME/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${reference} \
        -I ${bam} \
        -known ${known_indels_1} \
        -known ${known_indels_2} \
        -o ${sample_id}.indel_target_intervals.list
    """
}