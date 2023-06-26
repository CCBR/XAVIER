# exome-seek 🔬 [![tests](https://github.com/CCBR/exome-seek/workflows/tests/badge.svg)](https://github.com/CCBR/exome-seek/actions/workflows/main.yaml) [![docs](https://github.com/CCBR/exome-seek/workflows/docs/badge.svg)](https://github.com/CCBR/exome-seek/actions/workflows/docs.yml) [![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_wes_base)](https://hub.docker.com/r/nciccbr/ccbr_wes_base) [![GitHub issues](https://img.shields.io/github/issues/CCBR/exome-seek?color=brightgreen)](https://github.com/CCBR/exome-seek/issues)  [![GitHub license](https://img.shields.io/github/license/CCBR/exome-seek)](https://github.com/CCBR/exome-seek/blob/main/LICENSE) 

> **_GATK4 Whole Exome-sequencing Pipeline_**. This is the home of the pipeline, exome-seek. Its long-term goals: to accurately call germline and somatic variants, to infer CNVs, and to boldly annotate variants like no pipeline before!

## Overview
Welcome to exome-seek! Before getting started, we highly recommend reading through [exome-seek's documentation](https://CCBR.github.io/exome-seek).

The **`./exome-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>exome-seek <b>run</b></code>](https://CCBR.github.io/exome-seek/usage/run/): Run the GATK4 WES pipeline with your input files.
 * [<code>exome-seek <b>unlock</b></code>](https://CCBR.github.io/exome-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>exome-seek <b>cache</b></code>](https://CCBR.github.io/exome-seek/usage/cache/): Cache remote resources locally, coming soon!

Exome-seek is a comprehensive whole exome-sequencing pipeline following the Broad's set of best practices. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ or BAM files and can be run locally on a compute instance, on-premise using a cluster, or on the cloud (feature coming soon!). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://CCBR.github.io/exome-seek/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/CCBR/exome-seek/issues).

## Dependencies
**Requires:** `singularity>=3.5`  `snakemake==6.X`

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step relies on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singaularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity are the only two dependencies.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/CCBR/exome-seek.git
# Change your working directory
cd exome-seek/
```

## Contribute 

This site is a living document, created for and by members like you. EXOME-seek is maintained by the members of CCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [repository](https://github.com/CCBR/exome-seek/pulls).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
