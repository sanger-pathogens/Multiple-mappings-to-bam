FROM python:2.7-slim

COPY dependencies/bam_filter.py /opt/bam_filter/bam_filter.py

# Set up the Python script alias
RUN echo "alias bam_filter.py='python2 /opt/bam_filter/bam_filter.py'" >> ~/.bashrc

CMD ["/bin/bash"]