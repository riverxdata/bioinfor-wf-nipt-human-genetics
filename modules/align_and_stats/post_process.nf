process POST_PROCESS {
    container 'quay.io/biocontainers/samtools:1.15.1--h6899075_1'
    label 'process_medium'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id) path("${sample_id}.sorted.rmdup.{bam,bai}")

    script:
    """
    # Sort BAM
    samtools sort -@ ${task.cpus} -O bam -o ${sample_id}.sorted.bam ${bam_file}
    echo "** bam sorted done **"

    # Remove duplicates
    samtools rmdup ${sample_id}.sorted.bam ${sample_id}.sorted.rmdup.bam
    echo "** rmdup done **"

    # Index final BAM
    samtools index ${sample_id}.sorted.rmdup.bam
    echo "** index done **"
    """
}
