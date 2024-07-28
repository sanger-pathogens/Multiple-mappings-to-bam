FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    python2

COPY dependencies/bcf_2_pseudosequence.py /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py

# Set up the Python script alias
RUN echo "alias bcf_2_pseudosequence.py='python2 /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py'" >> ~/.bashrc

CMD ["/bin/bash"]