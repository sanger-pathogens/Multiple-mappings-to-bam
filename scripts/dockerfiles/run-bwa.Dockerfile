FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    samtools

CMD ["/bin/bash"]