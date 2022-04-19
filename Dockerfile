FROM python:3.9-slim

# Install Circos
RUN apt-get update && \
    apt-get install -y circos liblist-moreutils-perl libgomp1

# Install MGCplotter
RUN pip install -U pip && \
    pip install mgcplotter --no-cache-dir

CMD ["/bin/bash"]
