FROM python:3.8

LABEL org.opencontainers.image.title="CellPhoneDB"
LABEL org.opencontainers.image.description="CellPhoneDB can be used to search for a particular ligand/receptor, or interrogate your own HUMAN single-cell transcriptomics data."
LABEL org.opencontainers.image.authors="CellPhoneDB Team"

SHELL ["/bin/sh", "-c"]

ENV DEBIAN_FRONTEND=noninteractive
ENV CELLPHONEDB_RELEASE_PATH=/opt/cellphonedb/releases
ENV R_OS_IDENTIFIER=debian-11
ENV R_VERSION=4.1.3
ENV LD_LIBRARY_PATH="/opt/R/${R_VERSION}/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-11-openjdk-amd64/lib/server"

# install OS packages
RUN apt-get update && \
    apt-get install -yq git wget gdebi-core build-essential software-properties-common

# # install R
RUN wget -q "https://cdn.rstudio.com/r/${R_OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb" && \
    gdebi -n r-${R_VERSION}_1_amd64.deb && \
    rm -rf gdebi r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript && \
    Rscript -e "install.packages(c('ggplot2','pheatmap'), repos='https://packagemanager.rstudio.com/all/latest')" 

# install CellphoneDB
COPY . /tmp/cellphonedb
RUN pip3 install /tmp/cellphonedb --no-cache-dir && \
    mkdir -p /opt/cellphonedb/releases && \
    cellphonedb database download

# cleanup
RUN apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/*

# smoke test
CMD ["cellphonedb","--help"]