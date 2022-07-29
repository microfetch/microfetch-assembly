#!/usr/bin/env bash

PATH=/opt/conda/bin:/opt/conda/condabin:/opt/conda/bin:/opt/conda/condabin:/opt/conda/envs/trimmomatic/bin:/opt/conda/envs/assembly/bin:/opt/conda/envs/quast/bin:/opt/conda/envs/confindr/bin:/opt/conda/envs/kat/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#cat /etc/environment

#if ! pgrep -x java >/dev/null;
#then
#  echo "! preg -x java >/dev/null is truthy"
#else
#  echo "! preg -x java >/dev/null is falsey"
#fi

#if ! [ -f /app/skip-pipeline ];
#then
#  echo "Marking pipeline as already run."
#  touch /app/skip-pipeline
#else
#  echo "Pipeline already run once for debugging."
#  exit 0
#fi

# Check whether no runs have happened or last run has finished
if ! [ -d "$WORK_DIR" ] || [ -f "$OUTPUT_DIR/api_interaction/complete.txt" ] || ! pgrep -x java >/dev/null;
then
  if ! [ -d "$WORK_DIR" ];
  then
    echo "New run, skipping cleaning step."
  else
      if ! pgrep -x java >/dev/null;
      then
        echo "Java not running (likely unclean exit), relaunching."
      fi
    echo "Cleaning: remove $WORK_DIR, $OUTPUT_DIR; run nextflow clean -f"
    rm -rf "$WORK_DIR"
    rm -rf "$OUTPUT_DIR"
    nextflow clean -f
  fi

  echo "Launching new assembly pipeline"
  nextflow run /app/main.nf \
  --adapter_file /app/adapters.fas \
  --qc_conditions /app/qc_conditions_nextera_relaxed_dev.yml \
  --fastq_pattern '*{R,_}{1,2}.f*q.gz' \
  --depth_cutoff 100 \
  --prescreen_file_size_check 12 \
  --confindr_db_path /confindr_database \
  --careful \
  --kmer_min_copy 3 \
  --skip_quast_summary \
  --skip_quast_multiqc \
  --skip_fastqc_multiqc \
  --output_dir "$OUTPUT_DIR" \
  --api_url "$API_URL" \
  -w "$WORK_DIR" \
  -qs 1000 \
  -profile test

#else
#  echo "Pipeline is already running, doing nothing."
fi