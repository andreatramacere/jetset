.. JetSeT documentation master file

Installation
============

Requirements
~~~~~~~~~~~~~~~~~~~~~~

The following python packages are required:

- python 2.7 or >=3.6 (python 3 is suggested, python 2 should work still fine but it is not supported anymore)
- setuptools
- scipy
- numpy
- astropy
- matplotlib
- swig
- future
- iminuit
- corner
- six
- emcee
- pyyaml



A C compiler is also necessary, plus the SWIG wrapper generator.

All the dependencies can be installed following the Anaconda method
(suggested) **OR** the pip method, as described in the following sections.


I suggest to use anaconda and python3 (https://www.anaconda.com/download/)

Download the JeTseT code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the source code from:

https://github.com/andreatramacere/jetset/archive/stable.tar.gz

Uncompress the  archive

.. code-block:: bash

   jetset-stable.tar.gz



cd to  the dir

.. code-block:: bash

    cd jetset-stable


Installation using Anaconda (suggested method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the requiremets:


.. code-block:: bash

    conda install --file requirements.txt




Install JetSeT running on the command line:

.. code-block:: bash

        python setup.py clean
        python setup.py install

**run all the examples outside of the installation dir**


Installation using PIP
~~~~~~~~~~~~~~~~~~~~~~


Install requirements running on the command line:

.. code-block:: bash

    pip install -r requirements.txt ``

**if pip fails to install swig you can try one of the following methods:**

    - on linux Ubuntu:

    .. code-block:: bash

        - sudo apt-get install python-dev
        - sudo apt-get install swig


    - on linux Debian:

    .. code-block:: bash

        - sudo aptitude install python-dev
        - sudo aptitude install swig


    - on linux Fedora:

    .. code-block:: bash

        - sudo yum install python-dev
        - sudo yum install swig


    - on mac:

    .. code-block:: bash

        - ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
        - brew install swig


Install JetSeT running on the command line:

.. code-block:: bash

        python setup.py clean
        python setup.py install

**run all the examples outside of the installation dir**