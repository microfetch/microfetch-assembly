def version_message(String version) {
    println(
        """
        ==============================================
              Assembly Pipeline version ${version}
        ==============================================
        """.stripIndent()
    )
}

def help_message() {
    println(
        """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --input_dir /path/to/reads --fastq_pattern '*{R,_}{1,2}*.fastq.gz' --output_dir /path/to/output \
                             --other_param_1 --other_param_2 -resume
        Mandatory arguments:
        --input_dir      Path to input directory containing the fastq files to be assembled
        --fastq_pattern  The regular expression that will match fastq files e.g '*{R,_}{1,2}*.fastq.gz'
        --output_dir     Path to output directory
        --adapter_file   The path to a fasta file containing adapter sequences to trim from reads
                         (You can use the one supplied using --adapter_file adapaters.fas)
        Optional arguments:
        --single_end Add if only single rather than paired reads are available
        --prescreen_genome_size_check Size in bp of the maximum estimated genome to assemble. Without this any size genome assembly will be attempted
        --prescreen_file_size_check Minumum size in Mb for the input fastq files. Without this any size of file will be attempted (this and prescreen_genome_size_check are mutually exclusive)
        --cutadapt Whether to use the extra trimminmg step to remove 3' adapater sequences using cutadapt
        --depth_cutoff The estimated depth to downsample each sample to. If not specified no downsampling will occur
        --careful Turn on the SPAdes careful option which improves assembly by mapping the reads back to the contigs
        --minimum_scaffold_length The minimum length of a scaffold to keep. Others will be filtered out. Default 500
        --minimum_scaffold_depth The minimum depth of coverage a scaffold must have to be kept. Others will be filtered out. Default 3
        --qc_conditions Path to a YAML file containing pass/warning/fail conditions used by QualiFyr (https://gitlab.com/cgps/qualifyr)
        --full_output Output pre_trimming fastqc reports, merged_fastqs and corrected_fastqs. These take up signficant disk space
        """.stripIndent()
    )
}

def pipeline_start_message(String version, Map params){
    log.info "======================================================================"
    log.info "                  GHRU assembly pipeline"
    log.info "======================================================================"
    log.info "Running version   : ${version}"
    log.info "Fastq inputs      : ${params.input_dir}/${params.fastq_pattern}"
    log.info ""
    log.info "-------------------------- Other parameters --------------------------"
    params.sort{ it.key }.each{ k, v ->
        if (v){
            log.info "${k}: ${v}"
        }
    }
    log.info "======================================================================"
    log.info "Outputs written to path '${params.output_dir}'"
    log.info "======================================================================"

    log.info ""
}
def complete_message(Map params, nextflow.script.WorkflowMetadata workflow, String version){
    // Display complete message
    log.info ""
    log.info "Ran the workflow: ${workflow.scriptName} ${version}"
    log.info "Command line    : ${workflow.commandLine}"
    log.info "Completed at    : ${workflow.complete}"
    log.info "Duration        : ${workflow.duration}"
    log.info "Success         : ${workflow.success}"
    log.info "Work directory  : ${workflow.workDir}"
    log.info "Exit status     : ${workflow.exitStatus}"
    log.info ""
}

def error_message(nextflow.script.WorkflowMetadata workflow){
    // Display error message
    log.info ""
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}

def upload_error_report(Map params, nextflow.script.WorkflowMetadata workflow, String version){
    // Place the error report on Digital Ocean Spaces, and report failure to the API
    print "Uploading error report..."
    fileName = "${params.output_dir}/api_interaction/error.txt"
    errorFile = new File(fileName)
    errorFile.append(workflow.errorReport)
    call = "python ${workflow.projectDir}/templates/api_interaction/callback_api_error.py ${params.output_dir} ${params.api_url} ${fileName}"
    call.execute()
    println " complete."
}
