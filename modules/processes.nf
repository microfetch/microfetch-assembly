process PRE_SCREEN_GENOME_SIZE_ESTIMATION {
    tag { sample_id }

    input:
    tuple(val(sample_id), file(reads))

    output:
    tuple(val(sample_id), file('mash_stats.out'))

    """
    kat hist --mer_len 21  --thread 1 --output_prefix ${sample_id} ${reads[0]} > /dev/null 2>&1 \
    && minima=`cat  ${sample_id}.dist_analysis.json | jq '.global_minima .freq' | tr -d '\\n'`
    mash sketch -o sketch_${sample_id}  -k 32 -m \$minima -r ${reads[0]}  2> mash_stats.out
    """
}
