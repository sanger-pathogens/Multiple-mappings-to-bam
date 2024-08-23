# FROM quay.io/ssd28/gsoc-experimental/samtools:1.3
FROM quay.io/ssd28/gsoc-experimental/samtools:1.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    smalt

CMD ["/bin/bash"]