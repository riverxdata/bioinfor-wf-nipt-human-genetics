// modules 1
include { MAP_SE } from './modules/align_and_stats/align'
// include { POST_PROCESS } from './modules/align_and_stats/post_process'
// include { REALIGNER } from './modules/align_and_stats/re_align'
// include { RECABLIRATOR } from "./modules/align_and_stats/recablirate"
// include { STATS_BAM } from "./modules/align_and_stats/stats_bam"
// include { BEDTOOLS_GENOMECOV } from "./modules/align_and_stats/stats_cov"
// include { INDEX_BAM } from "./modules/align_and_stats/index_bam"

workflow {
    // inputs
    samples_ch = Channel
        .fromPath("/home/giangnguyen/Documents/dev/bioinfor-wf-nipt/data/samples/fq.list")
        .splitCsv(sep: '\t', header: false)
        .map { it -> 
        def (fq, sample_id, cell_id, lane_id) = it
        tuple(sample_id, cell_id, lane_id, file(fq),)
        }

    genome_index_ch = Channel.fromPath('/home/giangnguyen/Documents/dev/bioinfor-wf-nipt/data/ref_h38/Homo_sapiens_assembly38.fasta')
        .map { fasta ->
            def index_files = file("${fasta}.*")
            tuple(fasta, index_files)
        }

    // ch_dbsnp
    // ch_1k_gold
    // ch_known_index

    // module 1: Align and statistics
    MAP_SE(
        samples_ch,
        genome_index_ch
    )

    // POST_PROCESS(
    //     MAP_SE.out.aligned_bam
    // )

    // REALIGNER(
    //     POST_PROCESS.out
    //     .combine( ch_1k_gold )
    //     .combine( ch_known_index )
    // )
    
    // INDEX_BAM(
    //     REALIGNER.out
    // )

    // RECABLIRATOR(
    //     INDEX_BAM.out
    //     .combine( ch_dbsnp )
    //     .combine( ch_1k_gold )
    //     .combine( ch_known_index )
    // )

    // module 2
}

