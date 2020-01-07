.. install file

Installation
============

Install the JetSeT from ANACONDA (suggested)
--------------------------------------------
I suggest to use anaconda and python3 (https://www.anaconda.com/download/)

- create a virtual environment (not necessary):

.. code-block:: bash

   conda create --name jetset python=3.7 ipython jupyter

- install the code

.. code-block:: bash

   conda activate jetset


.. code-block:: bash

    conda install -c andreatramacere jetset


Install the JetSeT from source
------------------------------

Requirements
^^^^^^^^^^^^
The following python packages are required:
 - python 2.7 or >=3.6 (python 3 is suggested, python 2 should work still fine)
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

All the dependencies can be installed following the Anaconda method **OR** the pip method, as described below.

Download the code
^^^^^^^^^^^^^^^^^

- Get the source code from: https://gitlab.com/andreatramacere/jetset/-/archive/stable/jetset-stable.tar.gz
- Uncompress the  archive  `jetset-stable.tar.gz`
- cd to  the dir `jetset-stable`

Installation from source using Anaconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:

  .. code-block:: bash

    conda install --file requirements.txt

-  run on the command line

   .. code-block:: bash

       python setup.py clean

       python setup.py install

**run all the examples outside of the installation dir**

Installation from source using PIP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:

  .. code-block:: bash

    pip install -r requirements.txt `


- Install JetSeT: run on the command line:

  .. code-block:: bash

        python setup.py clean

        python setup.py install

**run all the examples outside of the installation dir**