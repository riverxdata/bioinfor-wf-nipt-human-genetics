process RECALIBRATOR {
    label 'process_medium'
    tag "$sample_id"
    container 'docker.io/pegi3s/gatk-3:3.8-0'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_dict
    path dbsnp
    path known_indels_1
    path known_indels_2

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.realign.BQSR.bam")

    script:
    """
    REF=\$(find -L ./ -name "*.fasta") 
    java -Xmx${task.memory.toMega()}m -jar /opt/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -nct ${task.cpus} \
        -R \$REF \
        -I ${bam}\
        --knownSites ${dbsnp[0]} \
        --knownSites ${known_indels_1[0]} \
        --knownSites ${known_indels_2[0]} \
        -o ${sample_id}.recal_data.table

    java -Xmx${task.memory.toMega()}m -jar /opt/GenomeAnalysisTK.jar \
        -T PrintReads \
        -nct ${task.cpus} \
        -R \$REF \
        -I ${bam}\
        --BQSR ${sample_id}.recal_data.table \
        -o ${sample_id}.sorted.rmdup.realign.BQSR.bam
    """
}