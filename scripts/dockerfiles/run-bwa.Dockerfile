FROM samtools-1.3:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa 
CMD ["/bin/bash"]