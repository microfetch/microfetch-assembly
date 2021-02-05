nextflow.enable.dsl=2
include {help_message; version_message; complete_message; error_message; pipeline_start_message} from './modules/messages'
include {default_params; check_params } from './modules/params_parser'
include {help_or_version} from './modules/params_utilities'
<<<<<<< HEAD
include {find_genome_size} from './modules/process_utilities'
include {find_total_number_of_bases} from './modules/process_utilities'
=======
include {find_genome_size; find_total_number_of_bases} from './modules/process_utilities'
>>>>>>> nigeria_develop

version = '2.0.0'

//  print help or version if required
// setup default params
default_params = default_params()
// merge defaults with user params
merged_params = default_params + params
// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)
// starting pipeline
pipeline_start_message(version, final_params)

<<<<<<< HEAD
include {GENOME_SIZE_ESTIMATION; PRE_SCREEN_FASTQ_FILESIZE; WRITE_OUT_FILESIZE_CHECK; DETERMINE_MIN_READ_LENGTH; QC_PRE_TRIMMING; TRIMMING; CUTADAPT; QC_POST_TRIMMING; FASTQC_MULTIQC; SPECIES_IDENTIFICATION; READ_CORRECTION; CHECK_FOR_CONTAMINATION; COUNT_NUMBER_OF_BASES; MERGE_READS} from './modules/processes' addParams(final_params)
=======
include {GENOME_SIZE_ESTIMATION; PRE_SCREEN_FASTQ_FILESIZE; WRITE_OUT_FILESIZE_CHECK; DETERMINE_MIN_READ_LENGTH; QC_PRE_TRIMMING; TRIMMING; CUTADAPT; QC_POST_TRIMMING; FASTQC_MULTIQC; SPECIES_IDENTIFICATION; READ_CORRECTION; CHECK_FOR_CONTAMINATION; COUNT_NUMBER_OF_BASES; MERGE_READS; SPADES_ASSEMBLY; FILTER_SCAFFOLDS} from './modules/processes' addParams(final_params)
>>>>>>> nigeria_develop

include {PRESCREEN_GENOME_SIZE_WORKFLOW; PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW} from './modules/workflows' addParams(final_params)


workflow {
    if (final_params.single_read){
        sample_id_and_reads = Channel
        .fromPath("${final_params.input_dir}/${final_params.fastq_pattern}")
        .map{ file -> tuple (file.baseName.replaceAll(/\..+$/,''), file)}
        .ifEmpty { error "Cannot find any reads matching: ${final_params.input_dir}/${final_params.fastq_pattern}" }

    } else {
        sample_id_and_reads = Channel
        .fromFilePairs("${final_params.input_dir}/${final_params.fastq_pattern}")
        .ifEmpty { error "Cannot find any reads matching: ${final_params.input_dir}/${final_params.fastq_pattern}" }
    }

    GENOME_SIZE_ESTIMATION(sample_id_and_reads)
    genome_sizes = GENOME_SIZE_ESTIMATION.out.map { sample_id, path -> find_genome_size(sample_id, path.text) }
    // pre-screen check based on genome size
    if (final_params.prescreen_genome_size_check) {
        sample_id_and_reads = PRESCREEN_GENOME_SIZE_WORKFLOW(genome_sizes, sample_id_and_reads)
    }
    // pre screen check based on file size
    if (final_params.prescreen_file_size_check){
        sample_id_and_reads = PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW(sample_id_and_reads)
    }
    // Assess read length and make MIN LEN for trimmomatic 1/3 of this value
    DETERMINE_MIN_READ_LENGTH(sample_id_and_reads)
    QC_PRE_TRIMMING(sample_id_and_reads)
    
    // CUTADAPT and QC_Post_Cutadapt
    if (final_params.cutadapt){
        CUTADAPT(sample_id_and_reads, final_params.adapter_file)
        min_trim_length_and_reads = DETERMINE_MIN_READ_LENGTH.out.join(CUTADAPT.out)
         // Trimmming step
        TRIMMING(min_trim_length_and_reads, final_params.adapter_file)
    } else {
        min_trim_length_and_reads = DETERMINE_MIN_READ_LENGTH.out.join(sample_id_and_reads)
        // Trimmming step
        TRIMMING(min_trim_length_and_reads, final_params.adapter_file)
    }
    //QC_post_Trimming
    QC_POST_TRIMMING(TRIMMING.out)
    // Multi QC
    FASTQC_MULTIQC(QC_POST_TRIMMING.out.fastqc_directories.collect())
    // Species ID
    SPECIES_IDENTIFICATION(TRIMMING.out)
    
    genome_size_trimmed_fastq = TRIMMING.out.join(genome_sizes)

    //Read Correction Step
    READ_CORRECTION(genome_size_trimmed_fastq)

    // Check for contamination
    CHECK_FOR_CONTAMINATION(READ_CORRECTION.out)

    ///base_counting
    if (final_params.single_read) {
        min_read_length_and_fastqs = DETERMINE_MIN_READ_LENGTH.out.join(READ_CORRECTION.out)
    } else {
        COUNT_NUMBER_OF_BASES(READ_CORRECTION.out)
        base_counts = COUNT_NUMBER_OF_BASES.out.map { sample_id, file -> find_total_number_of_bases(sample_id, file.text) }
        corrected_fastqs_and_genome_size_and_base_count = READ_CORRECTION.out.join(genome_sizes).join(base_counts).map{ tuple -> [tuple[0], tuple[1], tuple[2], tuple[3]]}
        MERGE_READS(corrected_fastqs_and_genome_size_and_base_count)
        min_read_length_and_fastqs = DETERMINE_MIN_READ_LENGTH.out.join(MERGE_READS.out)
    }

    SPADES_ASSEMBLY(min_read_length_and_fastqs)

    FILTER_SCAFFOLDS(SPADES_ASSEMBLY.out)    

    // >>>>>>>>>> COLOMBIA QUAST, QUAST_SUMMARY AND QUAST_MULTIQC PROCESSES HERE



}
