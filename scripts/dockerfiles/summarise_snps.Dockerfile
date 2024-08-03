FROM python:2.7-slim

RUN pip install --no-cache-dir \
    biopython \
    fisher \
    guppy3 \
    scipy \
    optparse

COPY dependencies/summarise_snps.py /opt/summarise_snps/summarise_snps.py

CMD ["/bin/bash"]