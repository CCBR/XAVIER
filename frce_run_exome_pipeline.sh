#/bin/bash

runmode=$1
output_dir=$2

# sh run_exome_pipeline.sh init /data/sevillas2/WES_hg38

module load snakemake/7.3.8
module load singularity/3.9.7

if [[ $runmode == "init" ]]; then
    ./exome-seek run \
    --runmode init \
    --input .tests/*.R?.fastq.gz \
    --output $output_dir \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
elif [[ $runmode == "dryrun" ]]; then
    ./exome-seek run \
    --runmode dryrun \
    --input .tests/*.R?.fastq.gz \
    --output $output_dir \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
elif [[ $runmode == "run" ]]; then
    ./exome-seek run \
    --runmode run \
    --input .tests/*.R?.fastq.gz \
    --output $output_dir \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
else
    echo "Runmode must be: init, dryrun, or run"
fi
