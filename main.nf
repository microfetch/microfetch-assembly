nextflow.enable.dsl=2
include {help_message; version_message; complete_message; error_message; pipeline_start_message} from './modules/messages'
include {default_params; check_params } from './modules/params_parser'
include {help_or_version} from './modules/params_utilities'
include {find_genome_size} from './modules/process_utilities'

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

include {PRE_SCREEN_GENOME_SIZE_ESTIMATION; WRITE_OUT_EXCLUDED_GENOMES; PRE_SCREEN_FASTQ_FILESIZE; WRITE_OUT_FILESIZE_CHECK; DETERMINE_MIN_READ_LENGTH; QC_PRE_TRIMMING} from './modules/processes' addParams(final_params)
workflow {
    reads = Channel
    .fromFilePairs("${params.input_dir}/${params.fastq_pattern}")
    .ifEmpty { error "Cannot find any reads matching: ${params.input_dir}/${params.fastq_pattern}" }

    // pre-screen check based on genome size
    if (final_params.prescreen_genome_size_check) {
        PRE_SCREEN_GENOME_SIZE_ESTIMATION(reads)
        genome_sizes = PRE_SCREEN_GENOME_SIZE_ESTIMATION.out.map { pair_id, file -> find_genome_size(pair_id, file.text) }
        // excluded genomes
        excluded_genomes_based_on_size = genome_sizes.filter { it[1] > final_params.prescreen_genome_size_check }
        WRITE_OUT_EXCLUDED_GENOMES(excluded_genomes_based_on_size)

        included_genomes_based_on_size = genome_sizes.filter { it[1] <= final_params.prescreen_genome_size_check }
        reads = reads.join(included_genomes_based_on_size).map { items -> [items[0], items[1]] }
    }
    // pre screen check based on file size
    if (final_params.prescreen_file_size_check){
        PRE_SCREEN_FASTQ_FILESIZE(reads)
        // filter files based on size
        // included genomes
        included_genomes_based_on_file_size = PRE_SCREEN_FASTQ_FILESIZE.out.filter { it[1].toFloat() >= final_params.prescreen_file_size_check }
        // excluded genomes
        excluded_genomes_based_on_file_size = PRE_SCREEN_FASTQ_FILESIZE.out.filter { it[1].toFloat() < final_params.prescreen_file_size_check }

        WRITE_OUT_FILESIZE_CHECK(PRE_SCREEN_FASTQ_FILESIZE.out)

        reads = reads.join(included_genomes_based_on_file_size).map { items -> [items[0], items[1]] } 
    }
    // Assess read length and make MIN LEN for trimmomatic 1/3 of this value
    DETERMINE_MIN_READ_LENGTH(reads)
    min_trim_length_and_reads = DETERMINE_MIN_READ_LENGTH.out.join(reads)
    QC_PRE_TRIMMING(reads)
}
