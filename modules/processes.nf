process GENOME_SIZE_ESTIMATION {
    tag { sample_id }

    input:
    tuple(val(sample_id), path(reads))

    output:
    tuple(val(sample_id), path('mash_stats.out'))

    script:
    if (params.kmer_min_copy)
      """
        mash sketch -o sketch_${sample_id}  -k 32 -m ${params.kmer_min_copy} -r ${reads[0]}  2> mash_stats.out
      """
    else
      """
      kat hist --mer_len 21  --thread 1 --output_prefix ${sample_id} ${reads[0]} > /dev/null 2>&1 \
      && minima=`cat  ${sample_id}.dist_analysis.json | jq '.global_minima .freq' | tr -d '\\n'`
      mash sketch -o sketch_${sample_id}  -k 32 -m \$minima -r ${reads[0]}  2> mash_stats.out
      """

    stub:
      """
      echo "Estimated genome size: 4.79477e+06" > mash_stats.out
      """
}

process WRITE_OUT_EXCLUDED_GENOMES {
    tag { sample_id }

    publishDir "${params.output_dir}/estimated_size_of_excluded_genomes"
    input:
    tuple(val(sample_id), val(genome_size))

    output:
    path("${sample_id}.estimated_genome_size.txt")

    script:
    """
    echo ${genome_size} > ${sample_id}.estimated_genome_size.txt
    """
}

process PRE_SCREEN_FASTQ_FILESIZE {
    tag { sample_id }

    input:
    tuple(val(sample_id), path(file_pair))

    output:
    tuple(val(sample_id), stdout)

    script:
    """
    stat -Lc %s ${file_pair[0]} |  awk '{printf "%f", \$1/1000000}'
    """

    stub:
    """
    echo 11000000
    """
}

process WRITE_OUT_FILESIZE_CHECK {
    tag { sample_id }

    input:
    tuple(val(sample_id), val(file_size))

    output:
    tuple(val(sample_id), path("${sample_id}.file_size_check.tsv"))

    script:
    """
    echo "file\tsize\n${sample_id}\t${file_size}" > ${sample_id}.file_size_check.tsv
    """
}
process DETERMINE_MIN_READ_LENGTH {
    tag { sample_id }

    input:
    tuple(val(sample_id), path(file_pair))

    output:
    tuple(val(sample_id), stdout)

    script:
    """
    gzip -cd ${file_pair[0]} | head -n 400000 | printf "%.0f" \$(awk 'NR%4==2{sum+=length(\$0)}END{print sum/(NR/4)/3}')
    """
}

process QC_PRE_TRIMMING {
    tag { sample_id }

    if (params.full_output){
        publishDir "${params.output_dir}/fastqc/pre_trimming",
        mode: 'copy',
        pattern: "*.html"
    }

    input:
    tuple(val(sample_id), path(file_pair))

    output:
    path('*.html')

    script:
    if (params.single_end){
        """
        fastqc -java=/opt/conda/envs/assembly/bin/java ${file_pair[0]}
        """
    } else {
        """
        fastqc -java=/opt/conda/envs/assembly/bin/java ${file_pair[0]} ${file_pair[1]}
        """
    }

    stub:
	    """
	    touch ${file_pair[0]}_fastqc.html
	    touch ${file_pair[1]}_fastqc.html
	    """
}

// TRIMMING
process TRIMMING {
  memory '4 GB'
  tag { sample_id }

  input:
  tuple(val(sample_id), val(min_read_length), path(reads))
  path('adapter_file.fas')

  output:
  tuple(val(sample_id), path('trimmed_fastqs/*.f*q.gz'))

  script:
  if (params.single_end) {
    method = "SE"
    file_input_and_outputs = "${reads[0]} trimmed_fastqs/${reads[0]}"
  } else {
    method = "PE"
    file_input_and_outputs = "${reads[0]} ${reads[1]} trimmed_fastqs/${reads[0]} /dev/null trimmed_fastqs/${reads[1]} /dev/null "
  }
  """
  mkdir trimmed_fastqs
  trimmomatic ${method} -threads 1 -phred33 ${file_input_and_outputs} ILLUMINACLIP:adapter_file.fas:2:30:10 SLIDINGWINDOW:4:20 LEADING:25 TRAILING:25 MINLEN:${min_read_length}
  """

  stub:
    """
    mkdir trimmed_fastqs
    touch trimmed_fastqs/${reads[0]}
    touch trimmed_fastqs/${reads[1]}
    """
}

// Post-Trimming QC
process QC_POST_TRIMMING {
  tag { sample_id }

  publishDir "${params.output_dir}/fastqc/post_trimming",
    mode: 'copy',
    pattern: "*.html"

  input:
  tuple(val(sample_id), path(reads) )

  output:
  path('*.html')
  tuple(val(sample_id), path("${sample_id}_R*_fastqc.txt"), emit: qc_post_trimming_files)
  path("*_fastqc_data"), emit: fastqc_directories

  script:
  if (params.single_end) {
    r1_prefix = reads[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
    """
    fastqc -java=/opt/conda/envs/assembly/bin/java ${reads[0]} --extract
    # rename files
    mv ${r1_prefix}_fastqc/summary.txt ${sample_id}_R1_fastqc.txt

    # move files for fastqc
    mkdir ${r1_prefix}_fastqc_data
    mv ${r1_prefix}_fastqc/fastqc_data.txt ${r1_prefix}_fastqc_data
    """
  } else {
  r1_prefix = reads[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
  r2_prefix = reads[1].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
  """
  fastqc -java=/opt/conda/envs/assembly/bin/java ${reads[0]} ${reads[1]} --extract
  # rename files
  mv ${r1_prefix}_fastqc/summary.txt ${sample_id}_R1_fastqc.txt
  mv ${r2_prefix}_fastqc/summary.txt ${sample_id}_R2_fastqc.txt

  # move files for fastqc
  mkdir ${r1_prefix}_fastqc_data
  mkdir ${r2_prefix}_fastqc_data
  mv ${r1_prefix}_fastqc/fastqc_data.txt ${r1_prefix}_fastqc_data
  mv ${r2_prefix}_fastqc/fastqc_data.txt ${r2_prefix}_fastqc_data
  """
  }

  stub:
	  r1_prefix = reads[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
	  r2_prefix = reads[1].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
    """
    touch ${sample_id}_R1_fastqc.txt
    touch ${sample_id}_R2_fastqc.txt
	  mkdir ${r1_prefix}_fastqc_data
	  mkdir ${r2_prefix}_fastqc_data
	  touch ${r1_prefix}_fastqc_data
    touch ${r2_prefix}_fastqc_data
    touch ${r1_prefix}_fastqc.html
    touch ${r2_prefix}_fastqc.html
    """
}

// Cutadapt
process CUTADAPT {
  tag { sample_id }
  publishDir "${params.output_dir}/pruned_fastqs",
  mode:'copy',
  pattern: '*.f*q.gz'

  input:
  tuple(val(sample_id), path(reads) )
  path('adapter_file.fas')

  output:
  tuple(val(sample_id), path('pruned_fastqs/*.f*q.gz') )

  script:
  if (params.single_end) {
    file_input_and_outputs = "-o pruned_fastqs/${reads[0]} ${reads[0]}"
  } else {
    file_input_and_outputs = "-o pruned_fastqs/${reads[0]}  -p pruned_fastqs/${reads[1]} ${reads[0]} ${reads[1]}"
  }
  """
  mkdir pruned_fastqs
  cutadapt -m 50 -j 2 -a file:'adapter_file.fas' ${file_input_and_outputs}
  """


}

//FastQC MultiQC
process FASTQC_MULTIQC {
  tag { 'multiqc for fastqc' }
  memory { 4.GB * task.attempt }

  publishDir "${params.output_dir}/quality_reports",
    mode: 'copy',
    pattern: "multiqc_report.html",
    saveAs: { "fastqc_multiqc_report.html" }

  when:
  ! params.skip_fastqc_multiqc

  input:
  path(fastqc_directories)

  output:
  path("multiqc_report.html")

  script:
  """
  multiqc --interactive .
  """
}

// Species ID
process SPECIES_IDENTIFICATION {
  tag { sample_id }

  input:
  tuple(val(sample_id), path(reads))  // from trimmed_fastqs_for_species_id

  output:
  tuple(val(sample_id), path("species_investigation*.tsv")) // into bactinspector_output

  script:
  """
  bactinspector check_species -fq ${reads[0]}
  """

  stub:
    """
    touch species_investigation_dummy.tsv
    """
}
// Read Corection
process READ_CORRECTION {
  tag { sample_id }

  if (params.full_output){
    publishDir "${params.output_dir}",
      mode: "copy",
      pattern: "corrected_fastqs/*.fastq.gz"
  }

  input:
  tuple(val(sample_id), path(reads), val(genome_size))

  output:
  tuple(val(sample_id), path("corrected_fastqs/*.f*q.gz") )

  script:
  if (params.single_end) {
    reads_argument = "-r ${reads[0]}"
  } else {
    reads_argument = "-r ${reads[0]} -r ${reads[1]}"
  }
  """
  mkdir corrected_fastqs
  lighter -od corrected_fastqs ${reads_argument} -K 32 ${genome_size}  -maxcor 1 2> lighter.out
  for file in corrected_fastqs/*.cor.fq.gz
  do
    new_file=\${file%.cor.fq.gz}.fastq.gz
    mv \${file} \${new_file}
  done
  """

  stub:
    """
    mkdir corrected_fastqs
    touch corrected_fastqs/ERR586796_1.fastq.gz
    touch corrected_fastqs/ERR586796_2.fastq.gz
    """
}


process CHECK_FOR_CONTAMINATION {
  tag {sample_id}
  cpus 2

  publishDir "${params.output_dir}/confindr",
    mode: 'copy',
    saveAs: { file -> "${sample_id}_${file}"}

  input:
  tuple(val(sample_id), path(file_pair))
  path(confindr_db_path)

  output:
  tuple(val(sample_id), path('confindr_report.csv'))

  script:
  def db_path = confindr_db_path.toString() == 'default_confindr_database' ? '/confindr_database' : confindr_db_path
  if (file_pair[0] =~ /_R1/){ // files with _R1 and _R2
    """
    confindr.py -i . -o . -d ${db_path} -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m
    """
  } else { // files with _1 and _2
    """
    confindr.py -i . -o . -d ${db_path} -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m -fid _1 -rid _2
    """
  }

	stub:
		"""
		cp /app/test_output/confindr/ERR586796_confindr_report.csv confindr_report.csv
		"""
}

process COUNT_NUMBER_OF_BASES {
  tag {sample_id}

  input:
  tuple(val(sample_id), path(reads))

  output:
  tuple(val(sample_id), path('seqtk_fqchk.out') )

  script:
  """
  seqtk fqchk -q 25 ${reads[0]} > seqtk_fqchk.out
  """
}

//Downsample reads
process DOWNSAMPLE_READS{
  tag { sample_id }

  input:
  tuple(val(sample_id), path(reads), val(genome_size), val(base_count))

  output:
  tuple(val(sample_id), path("downsampled_fastqs/*.f*q.gz") )

  script:
  if (params.depth_cutoff  && base_count/genome_size > params.depth_cutoff.toInteger()){
    downsampling_factor = params.depth_cutoff.toInteger()/(base_count/genome_size)
  } else {
    downsampling_factor = ""
  }
  """
  DOWNSAMPLING_FACTOR=${downsampling_factor}
  mkdir downsampled_fastqs
  for READ_FILE in ${reads}
  do
    if [[ -z \$DOWNSAMPLING_FACTOR ]]
    then
      # if no downsampling then just move input fastq to output location
      mv \${READ_FILE} downsampled_fastqs/\${READ_FILE}
    else
      seqtk sample  -s 12345 \${READ_FILE} ${downsampling_factor} | gzip > downsampled_fastqs/\${READ_FILE}
    fi
  done
  """
}
// merge reads
process MERGE_READS{
  tag { sample_id }

  if (params.full_output){
    publishDir "${params.output_dir}",
      mode: 'copy',
      pattern: "merged_fastqs/*.fastq.gz"
  }

  input:
  tuple(val(sample_id), path(reads))

  output:
  tuple(val(sample_id), path("merged_fastqs/*.f*q.gz") )

  script:
  """
  flash -m 20 -M 100 -t 1 -d merged_fastqs -o ${sample_id} -z ${reads[0]} ${reads[1]}
  """

  stub:
    """
    mkdir merged_fastqs
    cp /app/test_output/merged_fastqs/ERR586796* merged_fastqs
    """
}

process SPADES_ASSEMBLY {
  // memory { 4.GB * task.attempt }
  cpus { task.attempt == 1 ? 1 : 2 }

  tag { sample_id }

  input:
  tuple(val(sample_id), path(reads), val(min_read_length)) // from min_read_length_and_raw_fastqs

  output:
  tuple env(SPADES_SUCCESS), val(sample_id), path("final_fasta_dir/*.fasta") // into scaffolds

  script:
  spades_memory = 4 * task.attempt
  if (min_read_length.toInteger() < 25 ) { // this is the read length divided by 3 see trimming step
    kmers = '21,33,43,53'
  }
  else if (min_read_length.toInteger() < 50) {
    kmers = '21,33,43,53,63,75'
  } else {
    kmers = '21,33,43,55,77,99'
  }

  if (params.careful) {
    careful = "--careful"
  } else {
    careful = ""
  }

  if (params.single_end) {
    reads_argument = "--s1 ${reads}"
  } else {
    reads_argument = "--pe1-1 ${reads[1]} --pe1-2 ${reads[2]} --pe1-m ${reads[0]}"
  }

  if (task.attempt <=5 ) {
      """
      spades.py ${reads_argument} --only-assembler ${careful} -o . --tmp-dir /tmp/${sample_id}_assembly -k ${kmers} --threads 1 --memory ${spades_memory} || >&2 echo "SPAdes failed"
      rm -rf /tmp/${sample_id}_assembly
      mkdir final_fasta_dir
      if [[ -f "scaffolds.fasta" ]]; then
        SPADES_SUCCESS="TRUE"
        mv scaffolds.fasta final_fasta_dir/
      elif [[ -f "contigs.fasta" ]]; then
        SPADES_SUCCESS="TRUE"
        mv contigs.fasta final_fasta_dir/
      else
        SPADES_SUCCESS="FALSE"
        >&2 echo "No contigs found"
        exit 42
      fi
      """
   } else {
     """
     SPADES_SUCCESS="FALSE"
     mkdir final_fasta_dir
     touch final_fasta_dir/empty.fasta
     """
   }

	stub:
    """
    SPADES_SUCCESS="TRUE"
    mkdir final_fasta_dir
    cp /app/test_output/assemblies/pass/ERR586796.fasta final_fasta_dir/ERR586796.fasta
    """
}

// filter scaffolds to remove small and low coverage contigs
process FILTER_SCAFFOLDS {
  tag { sample_id }

  input:
  tuple(val(sample_id), path(scaffold_file))

  output:
  tuple(val(sample_id), path("${sample_id}.fasta"), emit: scaffolds_for_single_analysis)
  path("${sample_id}.fasta"), emit: scaffolds_for_combined_analysis

  script:
  """
  contig-tools filter -l ${params.minimum_scaffold_length} -c ${params.minimum_scaffold_depth} -f ${scaffold_file}

  if [[ -f "scaffolds.filter_gt_${params.minimum_scaffold_length}bp_gt_${params.minimum_scaffold_depth}.0cov.fasta" ]]
  then
    ln -s scaffolds.filter_gt_${params.minimum_scaffold_length}bp_gt_${params.minimum_scaffold_depth}.0cov.fasta ${sample_id}.fasta
  elif [[ -f "contigs.filter_gt_${params.minimum_scaffold_length}bp_gt_${params.minimum_scaffold_depth}.0cov.fasta" ]]
  then
    ln -s contigs.filter_gt_${params.minimum_scaffold_length}bp_gt_${params.minimum_scaffold_depth}.0cov.fasta ${sample_id}.fasta
  fi
  """

}

// assess assembly with Quast
process QUAST {
  tag { sample_id }

  publishDir "${params.output_dir}/quast",
    mode: 'copy',
    pattern: "report.tsv",
    saveAs: { file -> "${sample_id}_quast_" + file.split('\\/')[-1] }

  input:
  tuple(val(sample_id), path(contig_file))

  output:
  path("${sample_id}"), emit: quast_dir
  tuple(val(sample_id), path("${sample_id}/report.tsv"), emit: quast_report)

  """
  quast.py ${contig_file} -o .
  mkdir ${sample_id}
  ln -s \$PWD/report.tsv ${sample_id}/report.tsv
  """
}


// assess assembly with Quast but in a single file
process QUAST_SUMMARY {
  tag { 'quast summary' }
  memory { 4.GB * task.attempt }

  publishDir "${params.output_dir}/quast",
    mode: 'copy',
    pattern: "*report.tsv",
    saveAs: { file -> "combined_${file}"}

  when:
  ! params.skip_quast_summary

  input:
  path(contig_files)

  output:
  path("*report.tsv") optional true

  """
  quast.py --no-plots --no-html ${contig_files} -o .
  """
}


// QUAST MultiQC
process QUAST_MULTIQC {
  tag { 'multiqc for quast' }
  memory { 4.GB * task.attempt }

  publishDir "${params.output_dir}/quality_reports",
    mode: 'copy',
    pattern: "multiqc_report.html",
    saveAs: { "quast_multiqc_report.html" }

  when:
  ! params.skip_quast_multiqc

  input:
  path(quast_files)

  output:
  path("multiqc_report.html")

  script:
  """
  multiqc --interactive .
  """
}


process QUALIFYR {
  tag { sample_id }

  publishDir "${params.output_dir}/assemblies/pass",
    mode: 'copy',
    pattern: 'assemblies/pass/*',
    saveAs: { file -> file.split('\\/')[-1] }

  publishDir "${params.output_dir}/assemblies/warning",
    mode: 'copy',
    pattern: 'assemblies/warning/*',
    saveAs: { file -> file.split('\\/')[-1] }

  publishDir "${params.output_dir}/assemblies/failure",
    mode: 'copy',
    pattern: 'assemblies/failure/*',
    saveAs: { file -> file.split('\\/')[-1] }

  input:
  path(qc_conditions_yml)
  tuple(val(sample_id), path(quast_report), path(confindr_report), path(fastqc_reports), path(scaffold_file), path(bactinspector_report), path(file_size_check_output))

  output:
  path('assemblies/**/*')
  path("${sample_id}.qualifyr.json"), emit: json_files

	script:
  """
  # extract min and max genome sizes from bactinspector output, min file size from
  # file_size_check output and replace place holder in conditions file
  MAX_GENOME_LENGTH=\$(cat ${bactinspector_report} | awk -F'\t' 'NR == 2 {print \$8}')
  # if no species match set to 0
  if [ -z \$MAX_GENOME_LENGTH ]; then MAX_GENOME_LENGTH=0; fi
  # add wobble
  MAX_GENOME_LENGTH=\$(echo \$MAX_GENOME_LENGTH  | awk '{printf("%d",  \$1 * 1.1)}')

  MIN_GENOME_LENGTH=\$(cat ${bactinspector_report} | awk -F'\t' 'NR == 2 {print \$9}')
  # if no species match set to 0
  if [ -z \$MIN_GENOME_LENGTH ]; then MIN_GENOME_LENGTH=0; fi
  # add wobble
  MIN_GENOME_LENGTH=\$(echo \$MIN_GENOME_LENGTH  | awk '{printf("%d",  \$1 * 0.9)}')

  MIN_FILE_SIZE=\$(cat ${file_size_check_output} | awk -F'\t' 'NR == 2 {print \$2}')

  sed -i "s/MAX_GENOME_LENGTH/\${MAX_GENOME_LENGTH}/" ${qc_conditions_yml}
  sed -i "s/MIN_GENOME_LENGTH/\${MIN_GENOME_LENGTH}/" ${qc_conditions_yml}
  sed -i "s/MIN_FILE_SIZE/\${MIN_FILE_SIZE}/" ${qc_conditions_yml}


  result=`qualifyr check -y ${qc_conditions_yml} -f ${fastqc_reports} -c ${confindr_report}  -q ${quast_report} -b ${bactinspector_report} -z ${file_size_check_output} -s ${sample_id} 2> ERR`
  return_code=\$?
  if [[ \$return_code -ne 0 ]]; then
    exit 1;
  else
    if [[ \$result == "PASS" ]]; then
      qc_level="pass"
    elif [[ \$result == "WARNING" ]]; then
      qc_level="warning"
    elif [[ \$result == "FAILURE" ]]; then
      qc_level="failure"
    fi
    mkdir -p assemblies/\${qc_level}
    mv ${scaffold_file} assemblies/\${qc_level}/

    if [[ \$result != "PASS" ]]; then
      mv ERR assemblies/\${qc_level}/${sample_id}_qc_result.tsv
    fi
  fi


  # make json file
  qualifyr check -y ${qc_conditions_yml} -f ${fastqc_reports} -c ${confindr_report}  -q ${quast_report} -b ${bactinspector_report} -z ${file_size_check_output} -s ${sample_id} -j -o .
  """

  stub:
    """
    mkdir assemblies
    mkdir assemblies/pass
    mkdir assemblies/failure
    mkdir assemblies/warning
    cp /app/test_output/assemblies/pass/ERR586796.fasta assemblies/pass/ERR586796.fasta
    touch ${sample_id}.qualifyr.json
    """
}


process QUALIFYR_FAILED_SAMPLE {
  tag { sample_id }
  input:
  tuple(val(sample_id), val(file_size))
  tuple(path(quast_template), path(failed_sample_conditions_template), path(bactinspector_template), path(confindr_template), path(fastqc_template), path(file_size_check_template))

  output:
  path("${sample_id}.qualifyr.json")

  script:

  """
  sed -i "s/FILE_SIZE/${file_size}/" ${file_size_check_template}
  sed -i "s/MIN_FILE_SIZE/${params.prescreen_file_size_check}/" ${failed_sample_conditions_template}

  # make json file
  qualifyr check -y ${failed_sample_conditions_template} -f ${fastqc_template} ${fastqc_template} -c ${confindr_template}  -q ${quast_template} -b ${bactinspector_template} -z ${file_size_check_template} -s ${sample_id} -j -o .
  """

  stub:
    """
    touch ${sample_id}.qualifyr.json
    """
}

process QUALIFYR_REPORT {
  tag { 'qualifyr report' }

  publishDir "${params.output_dir}/quality_reports",
    mode: 'copy',
    pattern: "qualifyr_report.*"

  input:
  path(json_files)
  val(version)

  output:
  path("qualifyr_report.*")

  script:
  workflow_command = workflow.commandLine.replaceAll('"', '\\\\"')
//  """
//  qualifyr report -i . -c 'quast.N50,quast.# contigs (>= 0 bp),quast.# contigs (>= 1000 bp),quast.Total length (>= 1000 bp),quast.GC (%),confindr.contam_status,bactinspector.species' -s "Analysis with GHRU Assembly Pipeline version ${version}<br><br>Command line:<br>${workflow_command}"
//  """
  """
  qualifyr report -i . -c 'quast.N50,quast.# contigs,quast.Total length,quast.GC (%),confindr.contam_status,bactinspector.species' -s "Analysis with GHRU Assembly Pipeline version ${version}<br><br>Command line:<br>${workflow_command}"
  """

  stub:
		"""
		echo "quast.N50.metric_value\tquast.N50.check_result\n20\tPASS" > qualifyr_report.dummy.tsv
		"""

}

  process WRITE_ASSEMBLY_TO_DIR {
    tag { "assemblies to output" }

    publishDir "${params.output_dir}",
      mode: "copy"

    input:
    path(scaffold_files)

    output:
    path("assemblies")

    script:
    """
    mkdir assemblies
    mv ${scaffold_files} assemblies/
    """

  }

process REPORT_IGNORED_IDS {
  tag "report ignored ids"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
  val ignored_ids

  output:
  path("ignored_ids.txt")

  script:
  ignored_ids_string = ignored_ids.join(" ")
  """
  for IGNORED_ID in $ignored_ids_string
  do
      echo \$IGNORED_ID >> ignored_ids.txt
  done
  """
}

process GET_API_INPUT {
	// Read inputs from an API
	// This is designed to interface with the microfetch-pipeline API.
	// This initial call is to /request_assembly_candidate/ and should
	// receive a JSON file with the ENA record summary.
	// If no records are awaiting assembly, this will return status code 204.

	// TODO: accept candidate respecting API's upload_url field
	tag "$api_url"
	errorStrategy 'retry'
  maxRetries 3
  publishDir "${params.output_dir}/api_interaction", mode: 'copy'

	input:
		val api_url

	output:
		path("api_response.json")

	script:
		template "api_interaction/get_api_input.py"

// 	stub:
// 		"""
// 		echo GET ${api_url}request_assembly_candidate/
// 		echo Stubbed based on templates/api_response.json
// 		cp /app/templates/api_response.json api_response.json
// 		"""
}

process DOWNLOAD_FASTQ {
	// Extract the fastq links from an API json response
	// Download and zip files from fastq links
	tag "$json_file.baseName"
  publishDir "${params.output_dir}/api_interaction", mode: 'symlink'

	input:
		path json_file

	output:
		path('api_input_files/*.fastq.gz')

	script:
		"""
		#! /usr/bin/env python

		# adapted from https://stackoverflow.com/a/11768443

		from json import loads
		import shutil
		import urllib.request as request
		from contextlib import closing
		import os

		with open("${json_file}", "r") as f:
			j = loads(f.read())

		os.mkdir("api_input_files")

		sources = j['fastq_ftp'].split(';')
		for source in sources:
			print(f"Downloading {source}")
			os.system(f"wget --tries=5 --wait=2 --output-document api_input_files/{os.path.basename(source)} ftp://{source}")

		"""

// 	stub:
// 		"""
// 		mkdir api_input_files
// 		touch api_input_files/ERR586796_1.fastq.gz
// 		touch api_input_files/ERR586796_2.fastq.gz
// 		"""
}

process FILTER_ASSEMBLY_RESULT {
	debug true
	// Check whether the assembled genome is suitable for inclusion in AMR Watch
	// This runs a python script which uses a different QC JSON file for
	// each organism.
  publishDir "${params.output_dir}/api_interaction", mode: "copy"

	input:
		path(qc_file)
		path(api_response)

	output:
		path("post_assembly_filter.tsv")

	script:
		template "api_interaction/filter_assembly_result.py"

	stub:
		"""
		cp ${qc_file} post_assembly_filter.tsv
		"""
}

process UPLOAD_TO_SPACES {
	tag "$assembled_genome"
	// Upload an assembled genome and its assembly report to Spaces.
	input:
		path assembled_genome
		path qualifyr_report

	output:
		stdout

	script:
		template "api_interaction/upload_to_spaces.py"

// 	stub:
// 		"""
// 		echo ftp://example.com/assembly/pipeline/demo/genome.fasta
// 		"""

}

process CALLBACK_API {
	debug true
	// Send assembly details back to the server.
	// This interfaces with the microfetch-pipeline API.

	// TODO: post_assembly_filters
	tag "$api_url"
	errorStrategy 'retry'
  maxRetries 3

	input:
		val api_url
		val api_response
		val spaces_url
		val qualifyr_report

	script:
		template "api_interaction/callback_api.py"

// 	stub:
// 		"""
// 		echo PUT $api_url record/sample_id/
// 		echo assembly_result: success
// 		echo assembled_genome_url: $spaces_url
// 		echo qualifyr_report: [$qualifyr_report]
// 		"""
}
