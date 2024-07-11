# XAVIER - e**X**ome **A**nalysis and **V**ariant explor**ER** 🔬

[![tests](https://github.com/CCBR/XAVIER/workflows/tests/badge.svg)](https://github.com/CCBR/XAVIER/actions/workflows/main.yaml)
[![docs](https://github.com/CCBR/XAVIER/workflows/docs/badge.svg)](https://github.com/CCBR/XAVIER/actions/workflows/docs.yml)
[![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_wes_base)](https://hub.docker.com/r/nciccbr/ccbr_wes_base)
[![GitHub Issues](https://img.shields.io/github/issues/CCBR/XAVIER?color=brightgreen)](https://github.com/CCBR/XAVIER/issues)
[![MIT license](https://img.shields.io/github/license/CCBR/XAVIER)](https://github.com/CCBR/XAVIER/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12727316.svg)](https://doi.org/10.5281/zenodo.12727316)

> **_*XAVIER - eXome Analysis and Variant explorER*_**. This is the home of the pipeline, XAVIER. Its long-term goals: to accurately call germline and somatic variants, to infer CNVs, and to boldly annotate variants like no pipeline before!

## Overview

Welcome to XAVIER! Before getting started, we highly recommend reading through [xavier's documentation](https://CCBR.github.io/XAVIER).

The **`xavier`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions:

- [<code>xavier <b>run</b></code>](https://CCBR.github.io/XAVIER/usage/run/): Run the XAVIER pipeline with your input files.
- [<code>xavier <b>unlock</b></code>](https://CCBR.github.io/XAVIER/usage/unlock/): Unlocks a previous runs output directory.
- [<code>xavier <b>cache</b></code>](https://CCBR.github.io/XAVIER/usage/cache/): Cache remote resources locally, coming soon!

XAVIER is a comprehensive whole exome-sequencing pipeline following the Broad's set of best practices. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ or BAM files and can be run locally on a compute instance, on-premise using a cluster, or on the cloud (feature coming soon!). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://CCBR.github.io/XAVIER/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/CCBR/XAVIER/issues).

## Dependencies

**Requires:** `singularity>=3.5` `snakemake==6.X`

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step relies on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singaularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity are the only two dependencies.

## Run XAVIER pipeline

### Biowulf

```bash
# XAVIER is configured to use different execution backends: local or slurm
# view the help page for more information
module load ccbrpipeliner
xavier run --help

# @slurm: uses slurm and singularity execution method
# The slurm MODE will submit jobs to the cluster.
# It is recommended running XAVIER in this mode.

# Please note that you can dry-run the command below
# by providing  --runmode dryrun flag

# Do not run this on the head node!
# Grab an interactive node
sinteractive --mem=110g --cpus-per-task=12 --gres=lscratch:200
module load ccbrpipeliner

# First, initialize the output directory
xavier run \
--input data/*.R?.fastq.gz \
--output /data/$USER/xavier_hg38 \
--genome hg38 \
--pairs pairs.txt \
--targets Targets_hg38.bed \
--mode slurm \
--runmode init

# Second, do a dry run to visualize outputs
xavier run \
--input data/*.R?.fastq.gz \
--output /data/$USER/xavier_hg38 \
--genome hg38 \
--pairs pairs.txt \
--targets Targets_hg38.bed \
--mode slurm \
--runmode dryrun

# Then do a complete run
xavier run \
--input data/*.R?.fastq.gz \
--output /data/$USER/xavier_hg38 \
--genome hg38 \
--pairs pairs.txt \
--targets Targets_hg38.bed \
--mode slurm \
--runmode run
```

### FRCE

```bash
# grab an interactive node
srun --export all --pty --x11 bash

# add xavier to path correctly
. /mnt/projects/CCBR-Pipelines/pipelines/guis/latest/bin/setup

# add SIF cache and TMPDIR path
SIFCACHE="/mnt/projects/CCBR-Pipelines/SIFs/XAVIER"
TMPDIR="/scratch/cluster_scratch/$USER"

# run xavier

# Initialize and then dryrun (or run)
xavier run \
--input data/*.R?.fastq.gz \
--output /data/$USER/xavier_hg38 \
--genome hg38 \
--sif-cache $SIFCACHE \
--tmp-dir $TMPDIR \
--pairs pairs.txt \
--targets Targets_hg38.bed \
--mode slurm \
--runmode init # run

```

## Contribute

This site is a living document, created for and by members like you. XAVIER is maintained by the members of CCBR and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [repository](https://github.com/CCBR/XAVIER/pulls).

## References

<sup>**1.** Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.** Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>
