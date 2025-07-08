#BUILD FROM EXISITNG IMAGE (Cuda 11.6.22 NGC IMAGE)
FROM nvcr.io/nvidia/cuda:11.6.2-devel-ubuntu20.04

ENV DEBIAN_FRONTEND noninteractive
ARG CUDA_ARCHITECTURES="70;75;80;86;90"
ENV CUDA_ARCHITECTURES=${CUDA_ARCHITECTURES}

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        apt-transport-https \
        ca-certificates \
        dbus \
        fontconfig \
        gnupg \
        libasound2 \
        libfreetype6 \
        libglib2.0-0 \
        libnss3 \
        libsqlite3-0 \
        libx11-xcb1 \
        libxcb-glx0 \
        libxcb-xkb1 \
        libxcomposite1 \
        libxcursor1 \
        libxdamage1 \
        libxi6 \
        libxml2 \
        libxrandr2 \
        libxrender1 \
        libxtst6 \
        libgl1-mesa-glx \
        libxkbfile-dev \
        openssh-client \
        wget \
        xcb \
        xkb-data && \
    apt-get clean

RUN cd /tmp && \
    wget https://developer.nvidia.com/downloads/assets/tools/secure/nsight-systems/2024_5/nsight-systems-2024.5.1_2024.5.1.113-1_amd64.deb && \
    apt-get install -y ./nsight-systems-2024.5.1_2024.5.1.113-1_amd64.deb && \
    rm -rf /tmp/*

#MISSING KEY PROGRAMS
RUN apt-get update \
    && apt-get install -y \
    make \
    wget \
#    build-essential \
#    libstdc++6 \
#    libfmt-dev \
#    libssl-dev \
#    libssh-dev \
#    openssh-client \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /Software/

########################################
########## INSTALL LM V2.5 #############
########################################
#INSTALL MINICONDA3
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_22.11.1-1-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py38_22.11.1-1-Linux-x86_64.sh -b -p \
    && rm -f Miniconda3-py38_22.11.1-1-Linux-x86_64.sh

#Build LM
ADD Lattice-Microbes /Software/Lattice_Microbes_v2.5/
ADD Lattice-Microbes/lm_precomp.yml /Software/Lattice_Microbes_v2.5/lm_precomp.yml

RUN conda env update --file /Software/Lattice_Microbes_v2.5/lm_precomp.yml --prune
###pre-reqs
#RUN conda install -c conda-forge -y python=3.7.3 \
#    && conda install -c conda-forge -y git==2.33.1 gcc_linux-64==7.3.0 gxx_linux-64==7.3.0 binutils_linux-64==2.31.1 libstdcxx-ng==11.2.0 libgcc-ng==11.2.0 \
#    && conda install -c anaconda -y hdf5==1.10.4 h5py==2.10.0 protobuf==3.13.0.1 \
#orig
#    && conda install -c conda-forge -y numpy==1.19.2 cython==0.29.4 cmake==3.14.0 xlrd==1.2.0 swig==4.0.2 jupyter==1.0.0 matplotlib==3.4.3 ipywidgets==7.6.5 tqdm==4.62.3 pillow==8.3.2 jinja2==3.0.2 scipy==1.5.2 pybind11==2.8.1 pandas==1.3.4 pytables==3.5.2 biopython==1.79 
##remove jnotebook related stuff
##    && conda install -c conda-forge -y numpy==1.19.2 cython==0.29.4 cmake==3.14.0 xlrd==1.2.0 swig==4.0.2 matplotlib==3.4.3 tqdm==4.62.3 pillow==8.3.2 jinja2==3.0.2 scipy==1ipy.5.2 pybind11==2.8.1 pandas==1.3.4 pytables==3.5.2 biopython==1.79 

###copy files and build
#ADD Lattice_Microbes_v2.5_7.19.2023 /Software/Lattice_Microbes_v2.5/

RUN cd /Software/Lattice_Microbes_v2.5 \
    && mkdir build \
    && cd build \
    && cmake ../src -D MPD_GLOBAL_T_MATRIX=True -D MPD_GLOBAL_R_MATRIX=True \
    && make \
    && make install


