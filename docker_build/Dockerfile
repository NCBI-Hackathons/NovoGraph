FROM ubuntu:14.04

MAINTAINER Nathan Dunn <nathandunn@lbl.gov>
ENV DEBIAN_FRONTEND noninteractive

# Make sure the en_US.UTF-8 locale exists, since we need it for tests
RUN locale-gen en_US en_US.UTF-8 && dpkg-reconfigure locales

# FROM https://github.com/adamnovak/vg

RUN apt-get -qq update && \ 
    apt-get --no-install-recommends -y install maven2 openjdk-7-jdk \ 
    python3 mafft samtools && \ 
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


