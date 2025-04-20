FROM ubuntu:20.04
# run non-interactive
ENV DEBIAN_FRONTEND=noninteractive
# dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    cmake \
    g++ \
    make \
    libpthread-stubs0-dev \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev

    # && rm -rf /var/lib/apt/lists/*
WORKDIR /opt

# clone the BaseVar repository
RUN git clone https://github.com/ShujiaHuang/basevar.git -b v1.2.0
# build
WORKDIR /opt/basevar
RUN cd htslib && autoreconf -i && ./configure && make
RUN cd bin && g++ -std=c++17 -O3 -fPIC ../src/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar