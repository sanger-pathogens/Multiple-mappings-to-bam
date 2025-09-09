FROM debian:buster-slim

RUN apt-get update && \
    apt-get install -y procps

CMD ["/bin/bash"]