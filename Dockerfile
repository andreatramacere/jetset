#FROM continuumio/miniconda3

FROM python:3.7-slim
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook


ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}


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
RUN sudo apt-get install swig
RUN pip install -e https://github.com/andreatramacere/jetset#egg=jetset
ADD notebooks/QuickStart.ipynb $HOME/notebooks/QuickStart.ipynb



USER ${NB_USER}
ENV PATH /opt/conda/envs/jetset-env/bin:$PATH
#CMD echo "source activate jetset-env" > ~/.bashrc
WORKDIR /home/jovyan/notebooks
