def find_genome_size(sample_id, mash_output) {
  m = mash_output =~ /Estimated genome size: (.+)/
  genome_size = Float.parseFloat(m[0][1]).toInteger()
  return [sample_id, genome_size]
}

def find_total_number_of_bases(sample_id, seqtk_fqchk_ouput){
  m = seqtk_fqchk_ouput =~ /ALL\s+(\d+)\s/
  total_bases = new BigDecimal(m[0][1]).toLong() * 2 // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [sample_id, total_bases]
}