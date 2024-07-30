# Stage 1: Build environment
FROM debian:bullseye-20240612-slim AS build

RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake gcc g++ make build-essential autoconf automake libtool pkg-config wget ca-certificates tar bzip2 libncurses5-dev zlib1g-dev\
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and extract samtools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2/download -O samtools-1.3.tar.bz2 && \
    tar -xvjf samtools-1.3.tar.bz2 && \
    rm samtools-1.3.tar.bz2 && \
    cd samtools-1.3 && \
    ./configure && \
    make all all-htslib && \
    make install install-htslib

# Stage 2: Runtime environment
FROM debian:bullseye-20240612-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends libncurses5-dev zlib1g-dev\
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy samtools from the build stage
COPY --from=build samtools-1.3/samtools st/samtools

# Set up the alias
ENV PATH="/st:${PATH}"

CMD ["/bin/bash"]
