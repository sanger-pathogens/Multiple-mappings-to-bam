FROM python:2.7-slim

RUN pip install --no-cache-dir \
    Bio \
    pysam \
    optparse

COPY dependencies/join_dna_files_with_indels.py /opt/join_dna_files_with_indels/join_dna_files_with_indels.py

CMD ["/bin/bash"]