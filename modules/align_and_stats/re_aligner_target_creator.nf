process RE_ALIGNER_TARGET_CREATOR {
    label 'process_medium'
    tag "$sample_id"
    container 'docker.io/pegi3s/gatk-3:3.8-0'

    input:
    tuple val(sample_id), path(bam)
    path ref
    path ref_dict
    path known_indels_1
    path known_indels_2

    output:
    path "${sample_id}.indel_target_intervals.list"

    script:
    """
    REF=\$(find -L ./ -name "*.fasta") 
    java -Xmx${task.memory.toMega()}m -jar /opt/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R \${REF} \
        -I ${bam[0]} \
        -known ${known_indels_1[0]} \
        -known ${known_indels_2[0]} \
        -o ${sample_id}.indel_target_intervals.list
    """
}