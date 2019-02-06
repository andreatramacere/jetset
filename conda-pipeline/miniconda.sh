#!/bin/bash
################################################################################
##  File:  miniconda-install.sh
##  Team:  CI-Platform
##  Desc:  Installs miniconda
################################################################################


# Install Miniconda
    curl -sL https://repo.continuum.io/miniconda/Miniconda2-4.5.12-MacOSX-x86_64.sh -o miniconda.sh \
    && chmod +x miniconda.sh \
    && ./miniconda.sh -b -p /usr/share/miniconda \
    && rm miniconda.sh

CONDA=/usr/share/miniconda
echo "CONDA=$CONDA" | tee -a /etc/environment

ln -s $CONDA/bin/conda /usr/bin/conda

