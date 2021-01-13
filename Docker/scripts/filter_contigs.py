#!/usr/bin/env python
import sys
from Bio import SeqIO
import re
import argparse

cov_pattern = re.compile("cov_([0-9.]+)")

def filtered_contigs_generator(fasta_file_handle, min_length, min_coverage):
  total_contigs = 0
  contigs_kept = 0
  for contig in SeqIO.parse(fasta_file_handle, 'fasta'):
    total_contigs = total_contigs + 1
    result = cov_pattern.search(contig.name)
    if result:
      if float(result.group(1)) >= min_coverage and len(contig) >= min_length:
        contigs_kept = contigs_kept + 1
        yield contig
  print("Starting contigs: {0}\nContigs kept: {1}".format(total_contigs, contigs_kept))

def filter_contigs(fasta_file_path, min_contig_length, min_contig_coverage):
  with open(fasta_file_path.replace('fa', 'filter_gt_{0}bp_gt_{1}cov.fa'.format(min_contig_length, min_contig_coverage)), 'w') as filtered_fasta:
    with open(fasta_file_path) as input_fasta:
      SeqIO.write(filtered_contigs_generator(input_fasta, min_contig_length, min_contig_coverage), filtered_fasta, 'fasta')

def parse_arguments():
    description = """
    A function to filter contigs after assembly with SPAdes
    """
    # parse all arguments
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter, )
    parser.add_argument('-f', '--contigs_fasta_file', help='path to SPAdes contig fasta file', required = True, type = str)
    parser.add_argument('-l', '--minimum_contig_length', help='minimum length of a contig to keep', default = 500, type = int)
    parser.add_argument('-c', '--minimum_contig_coverage', help='minimum coverage of a contig to keep', default = 2.0, type = float)

    options = parser.parse_args()
    return options


if __name__ == "__main__":
    options = parse_arguments()
    filter_contigs(options.contigs_fasta_file, options.minimum_contig_length, options.minimum_contig_coverage)