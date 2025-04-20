process INDEL_RE_ALIGNER {
    label 'process_medium'
    tag "$sample_id"
    container 'docker.io/pegi3s/gatk-3:3.8-0'

    input:
    tuple val(sample_id), path(bam), path(target_interval)
    path ref
    path ref_dict
    path known_indels_1
    path known_indels_2

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.bam")

    script:
    """
    REF=\$(find -L ./ -name "*.fasta") 
    java -Xmx${task.memory.toMega()}m -jar /opt/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R \$REF \
        -I ${bam[0]} \
        -known ${known_indels_1[0]} \
        -known ${known_indels_2[0]} \
        --targetIntervals ${target_interval} \
        -o ${sample_id}.sorted.rmdup.realign.bam
    """
}