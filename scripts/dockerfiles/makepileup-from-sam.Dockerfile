FROM debian:bullseye-20240612-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bcftools \
    samtools \
    ca-certificates\
    openjdk-17-jre-headless

RUN mkdir /opt/picard
RUN wget -P /opt/picard https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
RUN echo "alias picard='java -jar /opt/picard/picard.jar'" >> ~/.bashrc

# Pending: Install GATK

CMD ["/bin/bash"]