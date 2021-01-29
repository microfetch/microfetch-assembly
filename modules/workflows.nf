include {WRITE_OUT_EXCLUDED_GENOMES; PRE_SCREEN_FASTQ_FILESIZE; WRITE_OUT_FILESIZE_CHECK} from './processes'

workflow PRESCREEN_GENOME_SIZE_WORKFLOW {
    take:
        genome_sizes
        sample_id_and_reads
    main:
        // excluded genomes
        excluded_genomes_based_on_size = genome_sizes.filter { it[1] > params.prescreen_genome_size_check }
        WRITE_OUT_EXCLUDED_GENOMES(excluded_genomes_based_on_size)

        included_genomes_based_on_size = genome_sizes.filter { it[1] <= params.prescreen_genome_size_check }
        included_sample_id_and_reads = sample_id_and_reads
                .join(included_genomes_based_on_size)
                .map { items -> [items[0], items[1]] }
    emit:
        included_sample_id_and_reads
}

workflow PRE_SCREEN_FASTQ_FILESIZE_WORKFLOW {
    take:
        sample_id_and_reads
    main:
        PRE_SCREEN_FASTQ_FILESIZE(sample_id_and_reads)
        // filter files based on size
        // included genomes
        included_genomes_based_on_file_size = PRE_SCREEN_FASTQ_FILESIZE.out.filter { it[1].toFloat() >= params.prescreen_file_size_check }
        // excluded genomes
        excluded_genomes_based_on_file_size = PRE_SCREEN_FASTQ_FILESIZE.out.filter { it[1].toFloat() < params.prescreen_file_size_check }

        WRITE_OUT_FILESIZE_CHECK(PRE_SCREEN_FASTQ_FILESIZE.out)

        included_sample_id_and_reads = sample_id_and_reads
            .join(included_genomes_based_on_file_size)
            .map { items -> [items[0], items[1]] }
    emit:
        included_sample_id_and_reads
}