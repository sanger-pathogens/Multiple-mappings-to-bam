FROM debian:bullseye-20240612-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    bcftools
CMD ["/bin/bash"]