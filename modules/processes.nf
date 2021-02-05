process GENOME_SIZE_ESTIMATION {
    tag { sample_id }

    input:
    tuple(val(sample_id), path(reads))

    output:
    tuple(val(sample_id), path('mash_stats.out'))

    script:
    """
    kat hist --mer_len 21  --thread 1 --output_prefix ${sample_id} ${reads[0]} > /dev/null 2>&1 \
    && minima=`cat  ${sample_id}.dist_analysis.json | jq '.global_minima .freq' | tr -d '\\n'`
    mash sketch -o sketch_${sample_id}  -k 32 -m \$minima -r ${reads[0]}  2> mash_stats.out
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
    if (params.single_read){
        """
        fastqc ${file_pair[0]}
        """
    } else {
        """
        fastqc ${file_pair[0]} ${file_pair[1]}
        """
    }
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
  if (params.single_read) {
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
  if (params.single_read) {
    r1_prefix = reads[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
    """
    fastqc ${reads[0]} --extract
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
  fastqc ${reads[0]} ${reads[1]} --extract
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
  if (params.single_read) {
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

  input:
  file(fastqc_directories) 

  output:
  file("multiqc_report.html")

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
  if (params.single_read) {
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
}


process CHECK_FOR_CONTAMINATION {
  tag {sample_id}
  cpus 2

  publishDir "${params.output_dir}/confindr",
    mode: 'copy',
    saveAs: { file -> "${sample_id}_${file}"}

  input:
  tuple(val(sample_id), path(file_pair))

  output:
  tuple(val(sample_id), path('confindr_report.csv'))

  script:
  if (file_pair[0] =~ /_R1/){ // files with _R1 and _R2
    """
    confindr.py -i . -o . -d /confindr_database -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m
    """
  } else { // files with _1 and _2
    """
    confindr.py -i . -o . -d  /confindr_database -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m -fid _1 -rid _2
    """  
  }

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
    downsampling_factor= null
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
}

process SPADES_ASSEMBLY {
  memory { 4.GB * task.attempt }
  
  tag { sample_id }

  input:
  tuple(val(sample_id), val(min_read_length), path(reads)) // from min_read_length_and_raw_fastqs
  
  output:
  tuple(val(sample_id), path("scaffolds.fasta")) // into scaffolds

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

  if (params.single_read) {
    reads_argument = "--s1 ${reads}"
  } else {
    reads_argument = "--pe1-1 ${reads[1]} --pe1-2 ${reads[2]} --pe1-m ${reads[0]}"
  }

  """
  spades.py ${reads_argument} --only-assembler ${careful} -o . --tmp-dir /tmp/${sample_id}_assembly -k ${kmers} --threads 1 --memory ${spades_memory}
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
  ln -s scaffolds.filter_gt_${params.minimum_scaffold_length}bp_gt_${params.minimum_scaffold_depth}.0cov.fasta ${sample_id}.fasta
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
  path("${sample_id}/report.tsv"), emit: quast_report 

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

  input:
  path(contig_files)

  output:
  path("*report.tsv") optional true

  """
  quast.py ${contig_files} -o .
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

  input:
  path(quast_files) 

  output:
  path("multiqc_report.html")

  script:
  """
  multiqc --interactive .
  """
}
