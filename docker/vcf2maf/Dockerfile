# Base image: VEP release 102 from Ensembl
# hub.docker.com/r/ensemblorg/ensembl-vep/tags
# VEP release 102 matches GENCODE v36 for GRCh38
# This GENCODE release matches the upcoming major
# relase to GDC GENCODE gene reference model.
# The GDC will be replacing the files in the 
# GDC data portal that were generated using 
# GENCODE v22 with files that were generated 
# using GENCODE v36. GDC will release the
# reprocessed data on 11/30/2021.
# @ubuntu/18.04
# @Dockerfile: 
# github.com/Ensembl/ensembl-vep/blob/release/102/docker/Dockerfile
FROM ensemblorg/ensembl-vep:release_102.0

LABEL maintainer=kuhnsa@nih.gov

# Create Container filesystem specific 
# working directory and opt directories
USER root
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system 
# packages required for your 
# project. If you need extra 
# system packages add them here.
# argparse.bash requires argparse
# python package.
# Installs python/3.6, which
# includes argparse by default
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        libbz2-dev \
        python3 \
        python3-pip \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install htslib/1.10
RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 \
        && tar -vxjf htslib-1.10.tar.bz2 \
        && rm htslib-1.10.tar.bz2 \
        && cd htslib-1.10 \
        && make \
        && cd /opt2

# Install samtools/1.10
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
        && tar -vxjf samtools-1.10.tar.bz2 \
        && rm samtools-1.10.tar.bz2 \
        && cd samtools-1.10 \
        && make \
        && cd /opt2

# Install bcftools/1.10
RUN wget https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2 \
        && tar -vxjf bcftools-1.10.tar.bz2 \
        && rm bcftools-1.10.tar.bz2 \
        && cd bcftools-1.10 \
        && make \
        && cd /opt2

# Get liftover v385 from UCSC
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/liftOver \
        && cd /opt2 \
        && chmod a+rx /opt2/liftOver

# Install vcf2maf/v1.6.21
RUN wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.21.tar.gz \
        && tar -xvf v1.6.21.tar.gz \
        && rm v1.6.21.tar.gz \
        && chmod a+rx /opt2/vcf2maf-1.6.21/*.pl


# Add Dockerfile and argparse.bash script
# and export environment variables
# and ensure proper permissions
ADD Dockerfile /opt2/Dockerfile
COPY argparse.bash /opt2
RUN chmod -R a+rX /opt2 \
    && chmod a+x /opt2/argparse.bash \
    && ln -sf /usr/bin/python3 /usr/bin/python
ENV PATH="$PATH:/opt2/bcftools-1.10:/opt2/samtools-1.10:/opt2/htslib-1.10:/opt2:/opt2/vcf2maf-1.6.21:/opt2/liftOver"

# Reset working directory
WORKDIR /data2