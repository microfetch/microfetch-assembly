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
        Mandatory arguments:
 
        Optional arguments:
        """.stripIndent()
    )
}

def pipeline_start_message(String version, Map params){
    log.info "======================================================================"
    log.info "                  GHRU assembly pipeline"
    log.info "======================================================================"
    log.info "Running version   : ${version}"
    log.info "Fastq inputs      : ${params.input_dir}/${params.fastq_pattern}"
    log.info "Adapter file      : ${params.adapter_file}"
    if (params.depth_cutoff){
        log.info "Depth cutoff      : ${params.depth_cutoff}"
    } else{
        log.info "Depth cutoff      : None"
    }

    log.info "======================================================================"
    log.info "Outputs written to path '${params.output_dir}'"
    log.info "======================================================================"
    log.info ""
}
def complete_message(Map params, nextflow.script.WorkflowMetadata workflow, String version){
    // Display complete message
    println ""
    println "Ran the workflow: ${workflow.scriptName} ${version}"
    println "Command line    : ${workflow.commandLine}"
    println "Completed at    : ${workflow.complete}"
    println "Duration        : ${workflow.duration}"
    println "Success         : ${workflow.success}"
    println "Work directory  : ${workflow.workDir}"
    println "Exit status     : ${workflow.exitStatus}"
    println ""
    println "Parameters"
    println "=========="
    params.sort{ it.key }.each{ k, v ->
        if (v){
            println "${k}: ${v}"
        }
    }
}

def error_message(nextflow.script.WorkflowMetadata workflow){
    // Display error message
    println ""
    println "Workflow execution stopped with the following message:"
    println "  " + workflow.errorMessage

}
