# Stage 1: Build Stage
FROM debian:bullseye-20240612-slim AS build

# Install dependencies for downloading and setting up GATK
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    procps \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and extract GATK
RUN mkdir /opt/gatk && \
    wget -P /opt/gatk https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip /opt/gatk/gatk-4.5.0.0.zip -d /opt/gatk && \
    rm /opt/gatk/gatk-4.5.0.0.zip

# Stage 2: Runtime Stage
# FROM quay.io/ssd28/gsoc-experimental/bcf_2_pseudosequence:0.0.1
FROM quay.io/ssd28/gsoc-experimental/bcf_2_pseudosequence:0.0.2

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    openjdk-17-jre-headless \
    python3 \
    python2 \
    procps \
    python-is-python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy GATK from the build stage
COPY --from=build /opt/gatk/gatk-4.5.0.0 /opt/gatk/gatk-4.5.0.0

# Set up the GATK alias
ENV PATH="/opt/gatk/gatk-4.5.0.0:${PATH}"

CMD ["/bin/bash"]
