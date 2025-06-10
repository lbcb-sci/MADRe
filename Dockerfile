FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Install dependencies and MADRe environment
RUN conda update -n base -c defaults conda && \
    conda create -y -n madre_env python=3.10 pip && \
    /bin/bash -c "source activate madre_env && \
    conda install -y -c bioconda -c conda-forge \
        scikit-learn minimap2 flye metamdbg hairsplitter seqkit && \
    pip install https://github.com/lbcb-sci/MADRe/archive/refs/tags/v1.0.0.tar.gz"

# Activate the environment automatically
SHELL ["/bin/bash", "-c"]
ENV PATH /opt/conda/envs/madre_env/bin:$PATH
ENV CONDA_DEFAULT_ENV madre_env

ENTRYPOINT ["madre"]
