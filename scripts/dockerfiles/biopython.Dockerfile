FROM python:2.7-slim

RUN apt-get update && \
    apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev procps

RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    numpy \
    biopython==1.68 

CMD ["/bin/bash"]