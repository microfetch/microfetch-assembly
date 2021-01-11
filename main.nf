nextflow.enable.dsl=2
include {help_message; version_message; complete_message; error_message; pipeline_start_message} from './lib/messages'
include {default_params; check_params } from './lib/params_parser'
include {help_or_version} from './lib/params_utilities'
include {PRE_SCREEN_GENOME_SIZE_ESTIMATION} from './lib/processes'

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

workflow {
    reads = Channel
    .fromFilePairs("${final_params.input_dir}/${final_params.fastq_pattern}")
    .ifEmpty { error "Cannot find any reads matching: ${final_params.input_dir}/${final_params.fastq_pattern}" }
    
    PRE_SCREEN_GENOME_SIZE_ESTIMATION(reads)
}
