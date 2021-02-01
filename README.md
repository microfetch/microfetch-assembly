### GHRU Assembly Pipeline for SE & PE Reads

Test command for paired end reads
```
nextflow run main.nf --input_dir small_test_input --output_dir test_output --fastq_pattern '*{R,_}{1,2}.fastq.gz' --adapter_file adapters.fas  --full_output --cutadapt -resume
```

Test command for single end reads
```
nextflow run main.nf --input_dir small_test_input --output_dir test_output --fastq_pattern '*{R,_}1.fastq.gz' --adapter_file adapters.fas  --full_output  --cutadapt --single_read -resume
```