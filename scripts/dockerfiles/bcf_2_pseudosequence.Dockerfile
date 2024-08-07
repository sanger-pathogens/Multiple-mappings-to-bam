FROM quay.io/ssd28/gsoc-experimental/samtools:1.3

RUN apt-get update && \
    apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev python2 pip
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools

RUN pip install --no-cache-dir \
    numpy \
    biopython==1.68

COPY dependencies/bcf_2_pseudosequence.py /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py

CMD ["/bin/bash"]