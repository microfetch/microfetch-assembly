process PRE_SCREEN_GENOME_SIZE_ESTIMATION {
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
  tuple(val(sample_id), path('trimmed_fastqs/*.f*q.gz') )

  script:
  if (params.single_read) {
    """
    mkdir trimmed_fastqs
    trimmomatic SE -threads 1 -phred33 ${reads[0]} trimmed_fastqs/${reads[0]} ILLUMINACLIP:adapter_file.fas:2:30:10 SLIDINGWINDOW:4:20 LEADING:25 TRAILING:25 MINLEN:${min_read_length}  

    """
  } else {
    """
    mkdir trimmed_fastqs
    trimmomatic PE -threads 1 -phred33 ${reads[0]} ${reads[1]} trimmed_fastqs/${reads[0]} /dev/null trimmed_fastqs/${reads[1]} /dev/null ILLUMINACLIP:adapter_file.fas:2:30:10 SLIDINGWINDOW:4:20 LEADING:25 TRAILING:25 MINLEN:${min_read_length}  
    """
  }
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
  tuple(val(sample_id), path("${sample_id}_R1_fastqc.txt"), path("${sample_id}_R2_fastqc.txt"), emit: qc_post_trimming_files)
  tuple(path("${r1_prefix}_fastqc_data"), path("${r2_prefix}_fastqc_data"), emit: fastqc_directories)

  script:
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

// Cutadapt
process CUTADAPT {
   memory '4 GB'
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

"""
mkdir pruned_fastqs
cutadapt -m 50 -j 2 -a file:"adapter_file.fas" -A file:"adapter_file.fas" -o pruned_fastqs/${reads[0]} -p pruned_fastqs/${reads[1]} ${reads[0]} ${reads[1]}

"""

}

// Post-Cutadapt QC
process QC_POST_CUTADAPT {
  tag { sample_id }

  publishDir "${params.output_dir}/fastqc/post_cutadapt",
    mode: 'copy',
    pattern: "*.html"
  
  input:
  tuple(val(sample_id), path(reads) )

  output:
  path('*.html')
  tuple(val(sample_id), path("${sample_id}_R1_fastqc.txt"), path("${sample_id}_R2_fastqc.txt"), emit: qc_post_trimming_files)
  tuple(path("${r1_prefix}_fastqc_data"), path("${r2_prefix}_fastqc_data"), emit: fastqc_directories)

  script:
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
