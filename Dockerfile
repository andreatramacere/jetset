#FROM continuumio/miniconda3

FROM python:3.7-slim
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook





USER root
#RUN conda create -n jetset-env python=3.7
#ENV PATH /opt/conda/envs/jetset-env/bin:$PATH
#RUN echo "source activate jetset-env" > ~/.bashrc

#RUN /bin/bash  source ~/.bashrc
#RUN /bin/bash  conda init bash
#RUN /bin/bash  conda activate jetset-env

ADD requirements_docker.txt /requirements_docker.txt
#RUN conda create -n jetset-env python=3.7 ipython notebook

#RUN conda install --yes -c astropy --file requirements_docker.txt
RUN pip install -r requirements_docker.txt
RUN apt-get update -y
RUN apt-get install -y swig
RUN apt-get install -y git
RUN apt-get install -y gcc
RUN pip install git+http://github.com/andreatramacere/jetset#egg=jetset
ADD notebooks/ $HOME/notebooks

# create user with a home directory
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER ${NB_USER}
WORKDIR ${HOME}/notebooks
