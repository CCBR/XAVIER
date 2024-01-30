# Base image
# Installs muTect/1.17 and java7 (jdk 1.7)
FROM ubuntu:14.04

LABEL maintainer=kuhnsa@nih.gov

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
# Strelka strictly requires python/2.7,
# it is not compatiable with any other
# version of python.
# Installs python/2.7.16
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      build-essential \
      bzip2 \
      git \
      locales \
      openjdk-7-jdk \
      software-properties-common \
      unzip \
      wget \
 && apt-get clean && apt-get purge \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install MuTect/1.1.7, from my mirror with a pre-compiled JAR file
# Requires java7 or jdk 1.7 for the pre-compiled JAR file (already satified)
# openjdk-7-jdk not a part of ubuntu >= 18.04, so installation using
# an older LTS version of Ubunutu is required
# Setting environment variables: MUTECT_JAR and MUTECT_JARPATH to
# mirror HPC module on Biowulf
RUN git clone https://github.com/skchronicles/mirror-muTect-1.1.7.git
ENV PATH="/opt2/mirror-muTect-1.1.7:$PATH"
ENV MUTECT_JAR="/opt2/mirror-muTect-1.1.7/muTect-1.1.7.jar"
ENV MUTECT_JARPATH="/opt2/mirror-muTect-1.1.7"

# Make java8 (jdk 1.8 the default version)
RUN ln -fs /usr/lib/jvm/java-7-openjdk-amd64/bin/java /usr/bin/java7

# Reset working directory
ADD Dockerfile /opt2/
RUN chmod -R a+rwX /opt2
WORKDIR /data2
