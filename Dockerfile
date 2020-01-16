FROM continuumio/miniconda3

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}


USER root
RUN conda create -n jetset-env python=3.7
ENV PATH /opt/conda/envs/jetset-env/bin:$PATH

RUN /bin/bash  source ~/.bashrc
RUN /bin/bash  conda init bash
RUN /bin/bash  conda activate jetset-env

ADD requirements_docker.txt /requirements_docker.txt
RUN conda create -n jetset-env python=3.7 ipython notebook

RUN conda install --yes -c astropy --file requirements_docker.txt
RUN conda install --yes -c andreatramacere -c astropy jetset
ADD notebooks/QuickStart.ipynb $HOME/notebooks



USER ${NB_USER}
ENV PATH /opt/conda/envs/jetset-env/bin:$PATH
CMD /bin/bash -c" source ~/.bashrc
CMD /bin/bash -c "conda init bash
CMD /bin/bash -c "conda activate jetset-env
WORKDIR /home/jovyan/notebooks
