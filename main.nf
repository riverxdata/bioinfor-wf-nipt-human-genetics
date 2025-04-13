// modules 1
include { ALIGN } from './modules/align_and_stats/align'
include { RE_ALIGNER_TARGET_CREATOR } from './modules/align_and_stats/re_aligner_target_creator'
include { INDEL_RE_ALIGNER } from './modules/align_and_stats/indel_re_aligner'
include { INDEX_BAM as INDEX_INDEL_REALIGN_BAM } from "./modules/align_and_stats/index_bam"
include { RECALIBRATOR } from "./modules/align_and_stats/recalibrator"
include { STATS_BAM } from "./modules/align_and_stats/stats_bam"
include { BEDTOOLS_GENOMECOV } from "./modules/align_and_stats/stats_cov"
include { BGZIP_AND_TABIX_BED } from "./modules/align_and_stats/bgzip_n_tabix_bed"
include { INDEX_BED } from "./modules/align_and_stats/index_bed"
// modules 2
include { GET_SNP_AND_AF } from "./modules/get_snp_n_af/get_snp_n_af"
// modules 3
include { PREPARE_REF_VCF } from "./modules/impute_genome_type/prepare_ref_vcf"
include { GET_IMPUTATION_CHUNK } from "./modules/impute_genome_type/get_imputation_chunk"
include { COMPUTE_GENOTYPE_LIKELIHOODS } from "./modules/impute_genome_type/compute_genotype_likelihoods"
include { MERGE_GENOTYPE_LIKELIHOODS } from "./modules/impute_genome_type/merge_genotype_likelihoods"
workflow {
    // inputs
    samples_ch = Channel
        .fromPath("/home/river/bioinfor-wf-nipt/data/samples/fq.list")
        .splitCsv(sep: '\t', header: false)
        .map { it -> 
        def (fq, sample_id, cell_id, lane_id) = it
        tuple(sample_id, cell_id, lane_id, file(fq),)
        }
    genome_index_ch = Channel
        .fromPath('/home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.fasta*', checkIfExists: true)
        .collect()
        .map { files ->
            def fasta_file = files.find { it.name.endsWith('.fasta') }
            def base_name = fasta_file.baseName
            tuple(base_name, files)
        }
    ref_ch = Channel
        .fromPath("/home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.{fasta,fasta.fai}",checkIfExists: true)
        .collect()
    
    ref_dict_ch = Channel        
        .fromPath("/home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.dict",checkIfExists: true)
        .collect()

    known_1k_gold_indel_ch = Channel
        .fromPath('/home/river/bioinfor-wf-nipt/data/gcs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz*', checkIfExists: true)
        .collect()
    
    known_assembly_indel_ch = Channel
        .fromPath('/home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.known_indels.vcf.gz*', checkIfExists: true)
        .collect()
    
    dbsnp_ch =  Channel
        .fromPath('/home/river/bioinfor-wf-nipt/data/gcs/Homo_sapiens_assembly38.dbsnp138.vcf*', checkIfExists: true)
        .collect()


    ref_panel_ch= Channel
    .fromPath('/home/river/bioinfor-wf-nipt/data/phasing/*.vcf.gz', checkIfExists: true)
    .map { vcf_file ->
        def match = (vcf_file.name =~ /(_chr[0-9XYM]+)\./)
        if (!match) exit 1, "Could not extract chromosome from ${vcf_file.name}"
        def chrom = match[0][1].replaceFirst('_', '') // get 'chr22' from '_chr22'
        def index_file = file("${vcf_file.toString()}.tbi")
        tuple(chrom, vcf_file, index_file)
    }

    // module 1: Align and statistic
    ALIGN(
        samples_ch,
        genome_index_ch
    )

    RE_ALIGNER_TARGET_CREATOR(
        ALIGN.out,
        ref_ch,
        ref_dict_ch,
        known_1k_gold_indel_ch,
        known_assembly_indel_ch
    )
    INDEL_RE_ALIGNER(
        ALIGN.out,
        ref_ch,
        ref_dict_ch,
        known_1k_gold_indel_ch,
        known_assembly_indel_ch,
        RE_ALIGNER_TARGET_CREATOR.out
    )

    INDEX_INDEL_REALIGN_BAM(
        INDEL_RE_ALIGNER.out
    )
    indel_realigned_bam_ch = INDEL_RE_ALIGNER.out.join(INDEX_INDEL_REALIGN_BAM.out)
    RECALIBRATOR(
        INDEL_RE_ALIGNER.out.join(INDEX_INDEL_REALIGN_BAM.out),
        ref_ch,
        ref_dict_ch,
        dbsnp_ch,
        known_1k_gold_indel_ch,
        known_assembly_indel_ch,
    )

    STATS_BAM(
        RECALIBRATOR.out
    )

    BEDTOOLS_GENOMECOV(
        RECALIBRATOR.out
    )

    BGZIP_AND_TABIX_BED(
        BEDTOOLS_GENOMECOV.out
    )
    // module 2
    GET_SNP_AND_AF(
        indel_realigned_bam_ch,
        ref_ch
    )
    // module 3
    // prepare reference panel and imputation chucks
    // it only runs once. We will devide it later
    PREPARE_REF_VCF(
        ref_panel_ch
    )
    
    GET_IMPUTATION_CHUNK(
        PREPARE_REF_VCF.out.vcf
    )

    COMPUTE_GENOTYPE_LIKELIHOODS(
        indel_realigned_bam_ch,
        PREPARE_REF_VCF.out.vcf,
        PREPARE_REF_VCF.out.tsv,
        ref_ch
    )
    grouped_vcfs_ch = COMPUTE_GENOTYPE_LIKELIHOODS.out
    .map { sample_id, chr, vcf, vcf_index ->
        tuple(chr, vcf, vcf_index)
    }
    .groupTuple().view()

    MERGE_GENOTYPE_LIKELIHOODS(
        grouped_vcfs_ch
    )
}