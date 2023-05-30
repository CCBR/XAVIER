#/bin/bash

runmode=$1

module load snakemake/7.19.1
module load singularity/3.10.5

if [[ $runmode == "init" ]]; then
    ./exome-seek run \
    --runmode init \
    --input .tests/*.R?.fastq.gz \
    --output /data/sevillas2/WES_hg38 \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
elif [[ $runmode == "dryrun" ]]; then
    ./exome-seek run \
    --runmode dryrun \
    --input .tests/*.R?.fastq.gz \
    --output /data/sevillas2/WES_hg38 \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
elif [[ $runmode == "run" ]]; then
    ./exome-seek run \
    --runmode run \
    --input .tests/*.R?.fastq.gz \
    --output /data/sevillas2/WES_hg38 \
    --genome hg38 \
    --targets .tests/Agilent_SSv7_allExons_hg38.bed
else
    echo "Runmode must be: init, dryrun, or run"
fi