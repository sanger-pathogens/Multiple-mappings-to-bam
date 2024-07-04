FROM debian:bullseye-20240612-slim

# Install dependencies, including Python, and remove unnecessary files in one layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bcftools \
    samtools \
    ca-certificates \
    openjdk-17-jre-headless \
    unzip \
    python3 \
    python2 \
    python-is-python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install GATK and clean up
RUN mkdir /opt/gatk && \
    wget -P /opt/gatk https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip /opt/gatk/gatk-4.5.0.0.zip -d /opt/gatk && \
    rm /opt/gatk/gatk-4.5.0.0.zip && \
    echo "alias gatk='/opt/gatk/gatk-4.5.0.0/gatk'" >> ~/.bashrc

COPY dependencies/bcf_2_pseudosequence.py /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py

RUN echo "alias bcf_2_pseudosequence.py='python2 /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py'" >> ~/.bashrc

CMD ["/bin/bash"]
