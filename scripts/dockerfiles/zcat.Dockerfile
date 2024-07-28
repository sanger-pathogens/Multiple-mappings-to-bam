FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    gzip
CMD ["/bin/bash"]