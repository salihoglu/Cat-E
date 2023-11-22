# Use an official R runtime as a parent image
#FROM rocker/r-ver:4.3.2

# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:latest

# System libraries
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    libglpk40
    
RUN sudo apt-get install default-jdk

RUN apt-get install -y xvfb x11vnc



RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev
    

# Copy Miniconda installer script
COPY shiny-app/Miniconda3-latest-Linux-x86_64.sh /opt/Miniconda3-latest-Linux-x86_64.sh

# Install Miniconda
RUN bash /opt/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm /opt/Miniconda3-latest-Linux-x86_64.sh

# Set conda in your PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Create a conda environment and activate it
RUN conda create --name cnapy-1.1.8 python=3.8

RUN pip install cobra

RUN pip install d3flux

RUN pip install playwright

RUN playwright install

RUN playwright install-deps


# Activate the conda environment
SHELL ["conda", "run", "-n", "cnapy-1.1.8", "/bin/bash", "-c"]

# Install Python packages within the conda environment
RUN conda install -n cnapy-1.1.8 -c conda-forge pandas numpy matplotlib scikit-learn

# Deactivate the conda environment
SHELL ["/bin/bash", "-c"]

# Install Bioconductor
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'


# Install gatom from GitHub using devtools
RUN R -e 'if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("gatom")'
    
RUN R -e 'if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("biomaRt")'

# Install Bioconductor packages
RUN R -e 'BiocManager::install(c("mwcsr", "limma", "edgeR", "SummarizedExperiment", "DiffCorr", "igraph","biomaRt","TCGAbiolinks","cBioPortalData","genefilter","GO.db","topGO","org.Hs.eg.db","clusterProfiler","GEOquery"))'

# Install other CRAN packages
RUN R -e 'install.packages(c("shiny","imsig","config","spsComps", "bs4Dash", "thematic", "waiter", "dplyr", "billboarder", "readxl", "shinyWidgets", "data.table", "visNetwork", "png", "plotly", "shinycssloaders", "ggplot2", "ggplotify", "DT", "RColorBrewer", "survival", "survminer", "reshape2", "maftools", "g3viz", "genefilter", "htmlwidgets", "shinyCyJS", "NGLVieweR", "httr", "colourpicker", "shinyBS", "rintrojs", "stringr", "ggvenn", "pathview", "heatmaply","RMySQL"))'

RUN mkdir /root/app

COPY /shiny-app /root/shiny_save


# Set up entry point script
COPY entrypoint.sh /usr/bin/entrypoint.sh
RUN chmod +x /usr/bin/entrypoint.sh

# Expose the Shiny app port
EXPOSE 3838

# Set the entry point to the script
ENTRYPOINT ["/usr/bin/entrypoint.sh"]

#EXPOSE 3838

# RUN dos2unix /usr/bin/shiny-server.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*
#CMD ["R", "-e", "shiny::runApp('/root/shiny_save', host='0.0.0.0', port=3838)"]

# RUN dos2unix /usr/bin/shiny-server.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*
#CMD ["R", "-e", "shiny::runApp('/shiny-app', host='0.0.0.0', port=3838)"]


#COPY /shiny-app ./shiny-app

# run app on container start
#CMD ["R", "-e", "shiny::runApp('/shiny-app')"]


