FROM python:2.7-slim

RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        procps \
        muscle &&
    apt-get clean

RUN pip install --upgrade pip && \
  pip install --no-cache-dir \
    numpy \
    biopython==1.68 \
    pysam==0.12.0.1

CMD ["/bin/bash"]