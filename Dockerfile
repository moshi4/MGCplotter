FROM python:3.9-slim

# Install MUMmer
RUN apt-get update && \
    apt-get install -y circos liblist-moreutils-perl

# Install GBKviz & Clear dependencies cache
RUN pip install -U pip && \
    pip install mgcplotter --no-cache-dir

CMD ["/bin/bash"]
