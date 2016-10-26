FROM ubuntu:16.04

# forked from https://github.com/edawson/vg_d
MAINTAINER Nathan Dunn <nathandunn@lbl.gov>
ENV DEBIAN_FRONTEND noninteractive

# Make sure the en_US.UTF-8 locale exists, since we need it for tests
RUN locale-gen en_US en_US.UTF-8 && dpkg-reconfigure locales

## VG START ##
# Install vg dependencies and clear the package index
RUN \
    echo "deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse" | tee -a /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
        build-essential \
        gcc-5-base \
        libgcc-5-dev \
        git \
        pkg-config \
        jq/trusty-backports

# Set up for make get-deps
RUN git clone --recursive https://github.com/vgteam/vg.git /app

RUN sed -i "s/sudo//g" /app/Makefile

   
# Move in all the other files
    
# Build vg
RUN cd /app && make get-deps && . ./source_me.sh && make -j8

# Make tests. We can't do it in parallel since it cleans up the test binary
RUN cd /app && . ./source_me.sh make test

ENV LD_LIBRARY_PATH=/app/lib
ENV LIBRARY_PATH /app/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /app/lib:$LD_LIBRARY_PATH
ENV LD_INCLUDE_PATH /app/include:$LD_INCLUDE_PATH
ENV C_INCLUDE_PATH /app/include:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH /app/include:$CPLUS_INCLUDE_PATH
ENV INCLUDE_PATH /app/include:$INCLUDE_PATH

#ENTRYPOINT ["/app/bin/vg"]
RUN cp /app/bin/vg /usr/bin


## VG END 

RUN apt-get -qq update && \ 
    apt-get --no-install-recommends -y install python3 mafft samtools parallel && \ 
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# RUN mkdir /app
COPY scripts /app/scripts/
COPY config /app/config/
COPY doc /app/doc/



