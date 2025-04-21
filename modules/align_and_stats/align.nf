process ALIGN {
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0'
    label 'process_medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), val(cell_id), val(lane_id), path(fq)
    tuple val(fasta), path(fasta_index)

    output:
    tuple val(sample_id), path("sorted.${sample_id}.{bam,bam.bai}"), emit: aligned_bam

    script:
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')

    bwa aln -e 10 -t ${task.cpus} -i 5 -q 0 \$INDEX ${fq} > ${sample_id}.sai

    bwa samse -r "@RG\\tID:${lane_id}\\tPL:COMPLETE\\tSM:${sample_id}" \\
        \$INDEX ${sample_id}.sai ${fq} | \\
        samtools sort -@ ${task.cpus} -O bam -o sorted.${sample_id}.bam -
    samtools index -@ ${task.cpus} sorted.${sample_id}.bam
    """
}
