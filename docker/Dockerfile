from ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y -q install \
    libgdal-dev \
    gdal-bin \
    python3-gdal \
    nano \
    vim \
    git \
    gcc \
    jags \
    python3-dev \
    python3-pip \
    python-setuptools \
    tzdata \
    software-properties-common \
    libffi-dev \ 
    libgeos-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libspatialindex-dev \
    fonts-dejavu \
    gfortran \
    wget && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install \
    azure-cli==2.14.0 \
    azure-storage-blob==12.5.0 \
    numpy==1.19.1 \
    pandas==1.1.0 \
    geopandas==0.8.1 \
    descartes==1.1.0 \
    shapely==1.7.0 \
    folium==0.11.0 \
    matplotlib==3.3.0 \
    rasterio==1.1.5 \
    dask[complete]==2.22.0 \
    rtree==0.9.4 \
    xarray==0.16.0 \
    seaborn==0.10.1 \
    jupyterlab==2.2.4 \
    progressbar==2.5 \
    progressbar2==3.53.1 \
    h5netcdf==0.8.1 \
    geopy==2.0.0 \
    zarr==2.4.0 \
    earthengine-api==0.1.231 \
    tqdm==4.48.2 \
    rasterstats==0.15.0 \
    gsutil==4.53 \
    tensorflow==2.3.0 \
    oauth2client==4.1.3 \ 
    google-auth-oauthlib==0.4.1 \
    netCDF4==1.5.4 \
    utm==0.6.0 \
    beautifulsoup4==4.10.0 \
    wget==3.2

#add vim bindings
RUN mkdir -p $(jupyter --data-dir)/nbextensions &&\
     cd $(jupyter --data-dir)/nbextensions &&\
     git clone https://github.com/lambdalisue/jupyter-vim-binding vim_binding &&\
     jupyter nbextension enable vim_binding/vim_binding

#add blobfuse dependencies
RUN wget https://packages.microsoft.com/config/ubuntu/18.04/packages-microsoft-prod.deb
RUN dpkg -i packages-microsoft-prod.deb
RUN apt-get update && apt-get install -y blobfuse

