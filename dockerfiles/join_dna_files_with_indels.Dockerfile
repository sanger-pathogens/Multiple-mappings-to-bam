FROM python:2.7-slim

RUN apt-get update && \
    apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev procps

RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    numpy \
    biopython==1.68 \
    pysam==0.12.0.1 

COPY dependencies/modules /opt/join_dna_files_with_indels/modules
COPY dependencies/join_dna_files_with_indels.py /opt/join_dna_files_with_indels/join_dna_files_with_indels.py

CMD ["/bin/bash"]