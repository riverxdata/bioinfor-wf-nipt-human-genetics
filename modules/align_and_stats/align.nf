process MAP_SE {
    container 'quay.io/biocontainers/bwa:0.7.17--pl5.22.0_2'
    label 'process_medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), val(cell_id), val(lane_id), path(fq)
    tuple val(fasta), path(index)

    output:
    tuple val(sample), path("${sample_id}.bam"), emit: aligned_bam
    script:
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')
    bwa aln -e 10 -t ${task.cpus} -i 5 -q 0 \$INDEX ${fq} > ${sample_id}.sai
    bwa samse -r "@RG\\tID:${lane_id}\\tPL:COMPLETE\\tSM:${sample_id}" \\
        \$INDEX ${sample_id}.sai ${fq} | \\
        samtools view -h -Sb - > ${sample_id}.bam
    """
}
