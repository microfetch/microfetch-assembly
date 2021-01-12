include {check_mandatory_parameter; check_optional_parameters} from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    params.help = false
    params.version = false
    params.input_dir = false
    params.output_dir = false
    params.fastq_pattern = false
    params.accession_number_file = false
    params.adapter_file = false
    params.depth_cutoff = false
    params.careful = false
    params.minimum_scaffold_length = false
    params.minimum_scaffold_depth = false
    params.confindr_db_path = false
    params.qc_conditions = false
    params.prescreen_genome_size_check = false
    params.prescreen_file_size_check = 5
    params.full_output = false
    return params
}

def check_params(Map params) { 
    final_params = params
    // check input dir
    final_params.input_dir = check_mandatory_parameter(params, 'input_dir') - ~/\/$/
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
    //  check a pattern has been specified
    if (params.input_dir){
        final_params.fastq_pattern = check_mandatory_parameter(params, 'fastq_pattern')
    }

    //  check an adapter_file has been specified
    final_params.adapter_file = check_mandatory_parameter(params, 'adapter_file')
    // // assign depth_cutoff from params
    // depth_cutoff = params.depth_cutoff

    // // set careful 
    // if (params.careful) {
    // careful = true
    // } else {
    // careful = false
    // }
    // // assign minimum scaffold length
    // if ( params.minimum_scaffold_length ) {
    // minimum_scaffold_length = params.minimum_scaffold_length
    // } else {
    // minimum_scaffold_length = 500
    // }
    // // assign minimum scaffold depth
    // if ( params.minimum_scaffold_depth ) {
    // minimum_scaffold_depth = params.minimum_scaffold_depth
    // } else {
    // minimum_scaffold_depth = 3
    // }

    // // set up confindr database path
    // if ( params.confindr_db_path ) {
    // confindr_db_path = params.confindr_db_path
    // } else {
    // // path in Docker image
    // confindr_db_path = "/home/bio/software_data/confindr_database"
    // }




    // // set full_output 
    // if (params.full_output) {
    // full_output = true
    // } else {
    // full_output = false
    // }
    return final_params
}

