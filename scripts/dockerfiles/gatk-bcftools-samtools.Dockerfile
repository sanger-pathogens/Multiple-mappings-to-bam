# Stage 1: Build Stage
FROM debian:bullseye-20240612-slim AS build

# Install dependencies for downloading and setting up GATK
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and extract GATK
RUN mkdir /opt/gatk && \
    wget -P /opt/gatk https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip /opt/gatk/gatk-4.5.0.0.zip -d /opt/gatk && \
    rm /opt/gatk/gatk-4.5.0.0.zip

# Stage 2: Runtime Stage
FROM samtools-1.3

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    openjdk-17-jre-headless \
    python3 \
    python2 \
    python-is-python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy GATK from the build stage
COPY --from=build /opt/gatk/gatk-4.5.0.0 /opt/gatk/gatk-4.5.0.0

# Set up the GATK alias
RUN echo "alias gatk='/opt/gatk/gatk-4.5.0.0/gatk'" >> ~/.bashrc

# Copy the Python script
COPY dependencies/bcf_2_pseudosequence.py /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py

# Set up the Python script alias
RUN echo "alias bcf_2_pseudosequence.py='python2 /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py'" >> ~/.bashrc

CMD ["/bin/bash"]
