FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    smalt \
    samtools

CMD ["/bin/bash"]