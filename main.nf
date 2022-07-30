nextflow.enable.dsl=2
// include non-process modules
include {help_message; version_message; complete_message; error_message; pipeline_start_message; upload_error_report} from './modules/messages'
include {default_params; check_params } from './modules/params_parser'
include {help_or_version} from './modules/params_utilities'
include {find_genome_size; find_total_number_of_bases; get_templates} from './modules/process_utilities'

version = '2.1.3'

// setup default params
default_params = default_params()
// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)
// starting pipeline
pipeline_start_message(version, final_params)

// include processes for pipelines
include {GENOME_SIZE_ESTIMATION; PRE_SCREEN_FASTQ_FILESIZE; WRITE_OUT_FILESIZE_CHECK; DETERMINE_MIN_READ_LENGTH; QC_PRE_TRIMMING; TRIMMING; CUTADAPT; QC_POST_TRIMMING; FASTQC_MULTIQC; SPECIES_IDENTIFICATION; READ_CORRECTION; CHECK_FOR_CONTAMINATION; COUNT_NUMBER_OF_BASES; DOWNSAMPLE_READS; MERGE_READS; SPADES_ASSEMBLY; FILTER_SCAFFOLDS; QUAST; QUAST_SUMMARY;  QUAST_MULTIQC; QUALIFYR; QUALIFYR_FAILED_SAMPLE; QUALIFYR_REPORT; WRITE_ASSEMBLY_TO_DIR; REPORT_IGNORED_IDS; GET_API_INPUT; DOWNLOAD_FASTQ; FILTER_ASSEMBLY_RESULT; UPLOAD_TO_SPACES; CALLBACK_API} from './modules/processes' addParams(final_params)
include {PRESCREEN_GENOME_SIZE_WORKFLOW; PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW} from './modules/workflows' addParams(final_params)

process DIE {
    """
    echo 'dying...'
    exit 1
    """
}

workflow {
		// set up input data
		if (final_params.api_url){
			// Read API input for fastq files, download them, and mock the fromFilePairs result
			api_json = GET_API_INPUT(final_params.api_url)
			DOWNLOAD_FASTQ(api_json)
			api_files = DOWNLOAD_FASTQ.out
			sample_id_and_reads = api_files
				.collect()
	      .map{ paths -> tuple(paths[0].baseName.replaceAll(/_[0-9+]\..+$/,''), tuple(paths[0], paths[1]))}
        .ifEmpty { exit 0, "No API input to process." }
		} else {
	    if (final_params.single_end){
	        sample_id_and_reads = Channel
	        .fromPath("${final_params.input_dir}/${final_params.fastq_pattern}")
	        .map{ file -> tuple (file.baseName.replaceAll(/\..+$/,''), file)}
	        .ifEmpty { error "Cannot find any reads matching: ${final_params.input_dir}/${final_params.fastq_pattern}" }
	    } else {
	        sample_id_and_reads = Channel
	        .fromFilePairs("${final_params.input_dir}/${final_params.fastq_pattern}")
	        .ifEmpty { error "Cannot find any reads matching: ${final_params.input_dir}/${final_params.fastq_pattern}" }
	    }
		}
		sample_id_and_reads.view()
    // Estimate genome sizes for pre screening (if specified) and for read correction
    GENOME_SIZE_ESTIMATION(sample_id_and_reads)
    genome_sizes = GENOME_SIZE_ESTIMATION.out.map { sample_id, path -> find_genome_size(sample_id, path.text) }
    // pre-screen check based on genome size
    if (final_params.prescreen_genome_size_check) {
        sample_id_and_reads = PRESCREEN_GENOME_SIZE_WORKFLOW(sample_id_and_reads)
    }
    // pre screen check based on file size
    if (final_params.prescreen_file_size_check){
        PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW(sample_id_and_reads)
        sample_id_and_reads = PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW.out[0]
        excluded_genomes_based_on_file_size = PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW.out[1]
        file_size_checks = PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW.out[2]
    }
    // Assess read length and make MIN LEN for trimmomatic 1/3 of this value
    DETERMINE_MIN_READ_LENGTH(sample_id_and_reads)
    // QC pre-trimming
    QC_PRE_TRIMMING(sample_id_and_reads)

    // Read Trimmimg
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
    //QC post-trimming
    QC_POST_TRIMMING(TRIMMING.out)
    // Multi QC
    FASTQC_MULTIQC(QC_POST_TRIMMING.out.fastqc_directories.collect())
    // Species ID
    SPECIES_IDENTIFICATION(TRIMMING.out)

    //Read Correction Step
    genome_size_trimmed_fastq = TRIMMING.out.join(genome_sizes)
    READ_CORRECTION(genome_size_trimmed_fastq)

    // Check for contamination
    CHECK_FOR_CONTAMINATION(READ_CORRECTION.out, final_params.confindr_db_path)

    corrected_reads = READ_CORRECTION.out
    // Downsample reads
    if (final_params.depth_cutoff){
        COUNT_NUMBER_OF_BASES(READ_CORRECTION.out)
        if (final_params.single_end) {
            base_counts = COUNT_NUMBER_OF_BASES.out.map { sample_id, file -> find_total_number_of_bases(sample_id, file.text, 1) }
        } else {
            base_counts = COUNT_NUMBER_OF_BASES.out.map { sample_id, file -> find_total_number_of_bases(sample_id, file.text, 2) }
        }
        corrected_fastqs_and_genome_size_and_base_count = READ_CORRECTION.out.join(genome_sizes).join(base_counts)
        DOWNSAMPLE_READS(corrected_fastqs_and_genome_size_and_base_count)
        corrected_reads = DOWNSAMPLE_READS.out
    }

    // Merge reads
    if (final_params.single_end) {
        min_read_length_and_fastqs = corrected_reads.join(DETERMINE_MIN_READ_LENGTH.out)
    } else {
        MERGE_READS(corrected_reads)
        min_read_length_and_fastqs = MERGE_READS.out.join(DETERMINE_MIN_READ_LENGTH.out)
    }
    // assemble reads
    SPADES_ASSEMBLY(min_read_length_and_fastqs)

    SPADES_ASSEMBLY.out
    .branch {
        success_sample_and_path  ->
        IDS_TO_IGNORE: success_sample_and_path[0] == "FALSE"
        IDS_TO_PROCESS: success_sample_and_path[0] == "TRUE"
    }
    .set { success_sample_and_paths_branch }

    // save ignored genomes in a different list
    success_sample_and_paths_branch.IDS_TO_IGNORE
        .map { it[1] }
        .set { ignored_ids }

    REPORT_IGNORED_IDS(ignored_ids.collect())

    // filter out small and low coverage contigs
    success_sample_and_paths_branch.IDS_TO_PROCESS
        .map { [it[1], it[2]] }
        .set { ids_to_filter }

    // filter out small and low coverage contigs
    FILTER_SCAFFOLDS(ids_to_filter)

    // run quast to assess quality of assemblies
    QUAST(FILTER_SCAFFOLDS.out.scaffolds_for_single_analysis)
    QUAST_SUMMARY(FILTER_SCAFFOLDS.out.scaffolds_for_combined_analysis.collect(sort: {a, b -> a.getBaseName() <=> b.getBaseName()}))
    QUAST_MULTIQC(QUAST.out.quast_dir.collect())

    // summarise quality
    if (final_params.qc_conditions){
		    quality_files = QUAST.out.quast_report
		      .join(CHECK_FOR_CONTAMINATION.out)
		      .join(QC_POST_TRIMMING.out.qc_post_trimming_files)
		      .join(FILTER_SCAFFOLDS.out.scaffolds_for_single_analysis)
		      .join(SPECIES_IDENTIFICATION.out)
		      .join(file_size_checks)
        QUALIFYR(final_params.qc_conditions, quality_files)
        QUALIFYR_FAILED_SAMPLE(excluded_genomes_based_on_file_size, get_templates())
        combined_qualifyr_json_files = QUALIFYR.out.json_files.concat(QUALIFYR_FAILED_SAMPLE.out).collect()
        // combined_qualifyr_json_files.view()
        QUALIFYR_REPORT(combined_qualifyr_json_files, version)
        if (final_params.api_url){
					UPLOAD_TO_SPACES(FILTER_SCAFFOLDS.out.scaffolds_for_combined_analysis, QUALIFYR_REPORT.out)
					FILTER_ASSEMBLY_RESULT(QUALIFYR_REPORT.out, GET_API_INPUT.out)
					CALLBACK_API(final_params.api_url, GET_API_INPUT.out, UPLOAD_TO_SPACES.out, QUALIFYR_REPORT.out, FILTER_ASSEMBLY_RESULT.out)
        }
    } else {
        scaffolds = FILTER_SCAFFOLDS.out.scaffolds_for_single_analysis.map{ tuple -> tuple[1]}.collect()
        WRITE_ASSEMBLY_TO_DIR(scaffolds)
				WRITE_ASSEMBLY_TO_DIR.out.view()
    }

}
workflow.onComplete {
    complete_message(final_params, workflow, version)
    "mkdir ${final_params.output_dir}/api_interaction".execute()
    "touch ${final_params.output_dir}/api_interaction/complete.txt".execute()
    println "Created ${final_params.output_dir}/api_interaction/complete.txt."
}

workflow.onError {
    error_message(workflow)
    "mkdir ${final_params.output_dir}/api_interaction".execute()
    "touch ${final_params.output_dir}/api_interaction/complete.txt".execute()
    println "Created ${final_params.output_dir}/api_interaction/complete.txt."
    if (final_params.api_url){
        upload_error_report(final_params, workflow, version)
    }
}