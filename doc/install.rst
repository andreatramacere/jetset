.. install file

Installation
============

.. important::
    Starting from version 1.1.0, python 2 is not supported anymore. Python >=3.6 is suggested, older python 3 versions (<3.6)  should work.


Install  JetSeT from Anaconda (suggested for Linux and Mac OSX)
-------------------------------------------------------------------





I suggest to use anaconda and python3 (https://www.anaconda.com/download/)

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

      conda create --name jetset python=3.7 ipython jupyter



  .. code-block:: bash

      conda activate jetset

- install the code

  .. code-block:: bash

    conda install -c andreatramacere -c astropy jetset



  if conda fails with dependencies you can try

  .. code-block:: bash

      conda install -c andreatramacere -c astropy jetset

  OR

  .. code-block:: bash

      conda install -c andreatramacere -c conda-forge jetset


- run the test

  .. code-block:: bash

      python -c 'from jetset.tests import test_functions; test_functions.test_short()'







Install the JetSeT from source
------------------------------


Download the code
^^^^^^^^^^^^^^^^^

- Get the source code from: https://github.com/andreatramacere/jetset/archive/stable.tar.gz
- Uncompress the  archive:  `jetset-stable.tar.gz`

- cd to  the dir source code dir

  .. code-block:: bash

   cd jetset-stable

Installation from source using Anaconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:

  .. code-block:: bash

    conda install --file requirements.txt

  if conda fails with dependencies you can try

  .. code-block:: bash

      conda install -c astropy --file requirements.txt

  OR

  .. code-block:: bash

      conda install -c conda-forge --file requirements.txt

-  run on the command line

   .. code-block:: bash

       python setup.py clean

       python setup.py install

- run the test (**run all the examples outside of the installation dir**)

   .. code-block:: bash

       python -c 'from jetset.tests import test_functions; test_functions.test_short()'





Installation from source using PIP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:

  .. code-block:: bash

    pip install -r requirements.txt `


- Install JetSeT: run on the command line:

  .. code-block:: bash

        python setup.py clean

        python setup.py install

- run the test  (**run all the examples outside of the installation dir**)

  .. code-block:: bash

       python -c 'from jetset.tests import test_functions; test_functions.test_short()'





Requirements
^^^^^^^^^^^^
The following python packages are required:
 - python >=3.6 (python >=3.6 is suggested, older python 3 versions should  work, python 2 is not supported any more from version>=1.1.0)
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

All the dependencies are installed following the Anaconda method **OR** the pip method, as described below.
