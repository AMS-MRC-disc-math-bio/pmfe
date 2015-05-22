# Build on a Debian "wheezy" foundation
FROM debian:jessie

MAINTAINER Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>

# Install the dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
            git \
            build-essential \
            gcc \
            make \
            libgmp-dev \
            libboost-filesystem-dev \
            libboost-program-options-dev \
            libboost-python-dev \
            libcgal-dev \
            ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Pull the pmfe repository and build the software
RUN cd /root && \
    git clone https://github.com/AMS-MRC-disc-math-bio/pmfe && \
    cd pmfe && \
    make -j
