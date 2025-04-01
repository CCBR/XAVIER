# Base image
FROM nciccbr/ccbr_ubuntu_base_20.04:v6

ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

LABEL maintainer <vishal.koparde@nih.gov, kuhnsa@nih.gov>

# Create Container filesystem specific
# working directory and opt directories
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2

# Set time zone to US east coast
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages required for your project
# If you need extra system packages add them here.
# python/3.8.0 and python/2.7.16 (strelka and manta)
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bc \
      build-essential \
      bzip2 \
      cmake \
      cpanminus \
      curl \
      g++ \
      gcc \
      gfortran \
      git \
      locales \
      make \
      openjdk-8-jdk \
      parallel \
      pigz \
      python2 \
      python3 \
      python3-pip \
      unzip \
      wget \
      less

# Make python3 the default interpreter
# and install Python Packages
#RUN ln -s /usr/bin/python3.8 /usr/bin/python
RUN pip3 install --upgrade pip \
	  && pip3 install argparse \
	  && pip3 install numpy \
      && pip3 install pysam \
	  && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Other misc dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      figlet \
      libatlas-base-dev \
      libblas-dev \
      libboost-dev \
      libbz2-dev \
      libexpat1-dev \
      libfreetype6-dev \
      libgd-dev \
      libgd-perl \
      libgs-dev \
      libgsl0-dev \
      libgsl-dev \
      libhtml-template-compiled-perl \
      libicu-dev \
      libjudy-dev \
      liblapack-dev \
      liblzma-dev \
      libmysqlclient-dev \
      libncurses-dev \
      libopenmpi-dev \
      libpng-dev \
      librtmp-dev \
      libssl-dev \
      libssl-dev \
      libtool \
      libxml2-dev \
      libxml-libxml-debugging-perl \
      libxml-opml-simplegen-perl \
      libxslt-dev \
      manpages-dev \
      zlib1g \
      zlib1g-dev \
      zlibc \
      openjdk-17-jdk \
      openjdk-11-jdk

# Install Java 17
# RUN apt-get update && apt-get install -y openjdk-17-jdk
# ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
# ENV PATH="$JAVA_HOME/bin:${PATH}"

# Common bioinformatics tools
# bwa/0.7.17  bowtie/1.2.3  bowtie2/2.3.5.1
# bedtools/2.27.1  bedops/2.4.37  samtools/1.10
# bcftools/1.10.2  vcftools/0.1.16  trimmomatic/0.39
# tabix/1.10.2
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bcftools \
      bedops \
      bedtools \
      bowtie \
      bowtie2 \
      bwa \
      samtools \
      tabix \
      trimmomatic \
      vcftools

# Install Sambamba/0.8.1 for CNTRL-Freec
# not available to apt-get on Ubuntu 20.04
RUN wget https://github.com/biod/sambamba/releases/download/v0.8.1/sambamba-0.8.1-linux-amd64-static.gz \
    && gzip -d /opt2/sambamba-0.8.1-linux-amd64-static.gz \
    && mv /opt2/sambamba-0.8.1-linux-amd64-static /opt2/sambamba \
    && chmod a+rx /opt2/sambamba

# Install GATK4 (GATK/4.6.0.0)
# Requires Java17
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip \
    && unzip /opt2/gatk-4.6.0.0.zip \
    && rm /opt2/gatk-4.6.0.0.zip \
    && /opt2/gatk-4.6.0.0/gatk --list
ENV PATH="/opt2/gatk-4.6.0.0:$PATH"

# Install last release of GATK3 (GATK/3.8-1)
# Only being used for the CombineVariants
# command that is not available in GATK4
# Available via env variable: $GATK_JAR
# Requires Java8 or 1.8
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && tar -xvjf /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && rm /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
ENV GATK_JAR="/opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

# Install dependencies needed to add a new repository over HTTPS
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      gnupg \
      dirmngr \
      ca-certificates \
      apt-transport-https \
      software-properties-common

# Install R (4.0) -- (edgeR needs 3.6) -- and R packages
# ggplot2  dplyr  plotly  htmlwidgets  tidyr and a few extras
# For more information, check out: https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev
RUN Rscript -e 'install.packages(c("argparse", "remotes", "knitr", "tidyverse", "dplyr", "plyr", "plotly", "ggplot2", "RColorBrewer", "htmlwidgets"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("shiny", "gridExtra", "flexdashboard", "rmarkdown", "crosstalk", "DT", "reshape2", "circlize", "viridis"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR", "ComplexHeatmap", "rtracklayer"))'

# Install Control-FREEC/v11.6 and additional dependencies
# Requires R, samtools, bedtools, sambamba (already satisfied)
# http://boevalab.inf.ethz.ch/FREEC/tutorial.html#install
RUN wget https://github.com/BoevaLab/FREEC/archive/refs/tags/v11.6.zip \
    && unzip /opt2/v11.6.zip \
    && rm /opt2/v11.6.zip \
    && cd /opt2/FREEC-11.6/src/ \
    && make
ENV PATH="/opt2/FREEC-11.6/src:$PATH"
WORKDIR /opt2

# Install Sequenza-Utils/3.0.0 and Sequenza
# Requires R, Python, SAMtools, tabix (already satisfied)
# https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html#getting-started
RUN pip3 install --upgrade pip \
	  && pip3 install sequenza-utils \
      && Rscript -e 'remotes::install_github("ShixiangWang/copynumber"); remotes::install_github("cran/sequenza")' \
      && Rscript -e 'remotes::install_version("iotools",version="0.3-2")'

# Install Somalier/v0.2.13
# download static binary
RUN wget https://github.com/brentp/somalier/releases/download/v0.2.13/somalier \
    && chmod a+rx /opt2/somalier

# Install VarScan/v2.4.6
# Works with java8
# and each wrapper script similar to HPC module
RUN wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.6.jar \
    && echo '#!/bin/bash' > /opt2/varscan \
    && echo 'java8 -jar /opt2/VarScan.v2.4.6.jar  "$@"' >> /opt2/varscan \
    && chmod a+rx /opt2/varscan

# Install snpEff/4.3t
# Works with java8
# Deleted bundled version of clinEff, not being used
# Setting env variable to jar file similar to HPC module
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip \
    && unzip /opt2/snpEff_v4_3t_core.zip \
    && rm /opt2/snpEff_v4_3t_core.zip \
    && rm -rf /opt2/clinEff/
ENV SNPEFF_JAR="/opt2/snpEff/snpEff.jar"

# Install samblaster/0.1.26
# Requires g++ and make (already satisfied)
RUN git clone https://github.com/GregoryFaust/samblaster.git \
    && cd /opt2/samblaster \
    && make \
    && chmod a+rX /opt2/samblaster/samblaster
ENV PATH="/opt2/samblaster:$PATH"
WORKDIR /opt2

# Install strelka/2.9.10, requires python2.7 (already satisfied)
# and set python2.7 as the default interpreter
RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && tar -xvjf /opt2/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && rm /opt2/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="/opt2/strelka-2.9.10.centos6_x86_64/bin:$PATH"

# Install manta/1.6.0
# Requires apt-get packages: bzip2 gcc g++ make python zlib1g-dev (already satisfied)
# Delete demo data (frees 275 MB)
RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && tar -xvjf /opt2/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm /opt2/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm -rf /opt2/manta-1.6.0.centos6_x86_64/share/demo/ \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="/opt2/manta-1.6.0.centos6_x86_64/bin:$PATH"

# Install VarDict/1.8.2
# Requires java8 (jdk 1.8 or later), perl and R (already satisfied)
RUN wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.2/VarDict-1.8.2.tar \
    && tar -xvf /opt2/VarDict-1.8.2.tar \
    && rm /opt2/VarDict-1.8.2.tar
ENV PATH="/opt2/VarDict-1.8.2/bin:$PATH"

# Add Dockerfile and argparse.bash script
# and export environment variables
# and creating a backwards compatible
# wrapper for Biowulf's trimmomatic/0.39
# and set java8 as default with links
# to alternative versions
ADD Dockerfile /opt2/base_gatk4_wes.dockerfile
COPY /argparse.bash /opt2
# Trimmomatic requires java11
ENV TRIMMOMATIC_JAR="/usr/share/java/trimmomatic-0.39.jar"
RUN echo '#!/bin/bash' > /opt2/trimmomatic \
    && echo 'java11 -jar "$TRIMMOMATIC_JAR" "$@"' >> /opt2/trimmomatic \
    && chmod +x /opt2/trimmomatic
RUN chmod -R a+rX /opt2 \
    && chmod a+x /opt2/argparse.bash \
    && ln -s /usr/lib/jvm/java-17-openjdk-amd64/bin/java /usr/bin/java17 \
    && ln -s /usr/lib/jvm/java-11-openjdk-amd64/bin/java /usr/bin/java11 \
    && ln -s /usr/lib/jvm/java-8-openjdk-amd64/bin/java /usr/bin/java8 \
    && ln -fs /usr/lib/jvm/java-17-openjdk-amd64/bin/java /usr/bin/java
ENV PATH="/opt2:$PATH"
ENV TRIMMOJAR="/usr/share/java/trimmomatic-0.39.jar"

WORKDIR /data2

# Clean-up step to reduce size
# and install GNU awk to calculate mean and standard
# deviation, ensures backward compatibility with
# biowulf installation of awk is a pointer to gawk,
# and install pandoc (>= 1.12.3 required for Rmarkdown)
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gawk \
    pandoc \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
