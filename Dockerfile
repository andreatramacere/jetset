FROM continuumio/miniconda3

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
ENV PATH /opt/conda/envs/env/bin:$PATH

USER root
RUN conda update
RUN conda create -n jetset-env python=3.7 ipython notebook
RUN echo "source activate jetset-env" > ~/.bashrc

RUN conda activate jetset-env
ADD requirements.txt /requirements_docker.txt
RUN conda install --yes -c astropy --file requirements_docker.txt
RUN conda install --yes -c andreatramacere -c astropy jetset
ADD notebooks/QuickStart.ipynb $HOME/notebooks



USER ${NB_USER}
CMD conda activate jetset-env
WORKDIR /home/jovyan/notebooks
