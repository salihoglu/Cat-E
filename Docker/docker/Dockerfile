# Use an official R runtime as a parent image
#FROM rocker/r-ver:4.3.2

# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:latest


# System libraries
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libnetcdf-dev \
    libcurl4-gnutls-dev \
    libhts-dev \ 
    libbz2-dev \ 
    liblzma-dev \ 
    libxml2-dev \
    python3 \
    python3-pip \
    libglpk40 \
    default-jdk \
    libmariadb-dev-compat \
    xvfb \
    x11vnc \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libmysqlclient21 \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    wget


# Install miniconda
#ENV CONDA_DIR /opt/conda
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod u+x ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Create and activate the conda environment
RUN conda create --name cnapy-1.1.8 -c Gurobi -c IBMDecisionOptimization -c conda-forge -c cnapy cnapy=1.1.8 --yes

RUN pip install cobra d3flux playwright rdkit
RUN playwright install --with-deps


# Activate the conda environment
SHELL ["conda", "run", "-n", "cnapy-1.1.8", "/bin/bash", "-c"]

# Install Python packages within the conda environment
RUN conda install -n cnapy-1.1.8 -c conda-forge pandas numpy matplotlib scikit-learn

# Deactivate the conda environment
SHELL ["/bin/bash", "-c"]

# Install Bioconductor
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'


# Install gatom from GitHub using devtools
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("gatom")'
    
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("biomaRt")'

# Install Bioconductor packages
RUN R -e 'BiocManager::install(c("mwcsr", "limma", "edgeR", "SummarizedExperiment","Rhtslib","DiffCorr","maftools", "igraph","cBioPortalData","genefilter","GO.db","pathview","topGO","org.Hs.eg.db","clusterProfiler","GEOquery"))'

# Install other CRAN packages
RUN R -e 'install.packages(c("shiny","imsig","reticulate","config","spsComps", "bs4Dash", "thematic", "waiter", "dplyr", "billboarder", "readxl", "shinyWidgets", "data.table", "visNetwork", "png", "plotly", "shinycssloaders", "ggplot2", "ggplotify", "DT", "RColorBrewer", "survival", "survminer", "reshape2", "g3viz",  "htmlwidgets", "shinyCyJS", "NGLVieweR", "httr", "colourpicker", "shinyBS", "rintrojs", "stringr", "ggvenn", "heatmaply","RMySQL"))'

# Check if packages are installed
RUN R -e 'if (!require("maftools", quietly = TRUE)) stop("Package maftools not found")'

RUN mkdir '/app'

COPY shiny-app/ app/


# Set up entry point script
COPY /docker/entrypoint.sh /usr/bin/entrypoint.sh
RUN chmod +x /usr/bin/entrypoint.sh

# Expose the Shiny app port
EXPOSE 3838

# Set the entry point to the script
ENTRYPOINT ["/usr/bin/entrypoint.sh"]
