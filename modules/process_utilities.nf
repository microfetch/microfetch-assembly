def find_genome_size(sample_id, mash_output) {
  m = mash_output =~ /Estimated genome size: (.+)/
  genome_size = Float.parseFloat(m[0][1]).toInteger()
  return [sample_id, genome_size]
}

def find_total_number_of_bases(sample_id, seqtk_fqchk_ouput, num_of_read_files){
  m = seqtk_fqchk_ouput =~ /ALL\s+(\d+)\s/
  total_bases = new BigDecimal(m[0][1]).toLong() * num_of_read_files // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [sample_id, total_bases]
}

def get_templates(){
  failed_sample_conditions_template = "${baseDir}/templates/failed_sample_qualifyr_files/conditions_for_failed_samples.yml"
  bactinspector_template = "${baseDir}/templates/failed_sample_qualifyr_files/bactinspector.tsv"
  confindr_template = "${baseDir}/templates/failed_sample_qualifyr_files/confindr.csv"
  fastqc_template = "${baseDir}/templates/failed_sample_qualifyr_files/fastqc.txt"
  file_size_check_template = "${baseDir}/templates/failed_sample_qualifyr_files/file_size_check.tsv"
  quast_template = "${baseDir}/templates/failed_sample_qualifyr_files/quast.txt"
  return([failed_sample_conditions_template, bactinspector_template, confindr_template, fastqc_template, file_size_check_template, quast_template])
}