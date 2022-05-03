.. _install_file:

Installation
============

.. important::
    Starting from version 1.1.0, python 2 is not supported anymore. Python >=3.6 is suggested, older python 3 versions (<3.6)  should work.


Windows 10 prerequisites
------------------------
Install the Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install-win10


Install  JetSeT from Anaconda
------------------------------------------------------------------------------
to get anaconda: https://www.anaconda.com/download/

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

      conda create --name jetset python=3.9 ipython jupyter

  .. code-block:: bash

      conda activate jetset

- install the code

  .. code-block:: bash

      conda install -c andreatramacere -c astropy -c conda-forge jetset



- run the test (optional)

  .. code-block:: bash

    pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser


Install  JetSeT from pip
------------------------------------------------------------------------------

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

    pip install virtualenv

    virtualenv -p python3.9 jetset

    source jetset/bin/activate

    pip install ipython jupyter



- install the code

  .. code-block:: bash

     pip install jetset



- run the test (optional)

  .. code-block:: bash

    pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser




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

      conda install --yes   swig">=3.0.0"

      conda install -c astropy -c conda-forge --file requirements.txt

.. important::
    if anaconda fails to install swig, you can try one of the following alternative :ref:`swig'

-  run on the command line

   .. code-block:: bash

       python setup.py clean

       python setup.py install

- run the test (optional, **run all the examples outside  the installation dir**)

  .. code-block:: bash

     cd ~/

     mkdir test_jetset

     cd test_jetset

     pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser






Installation from source using PIP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:

  .. code-block:: bash

    pip install swig>=3.0.0

    pip install -r requirements.txt

.. important::
    if pip fails to install swig, you can try one of the following alternative :ref:`swig'


- Install JetSeT: run on the command line:

  .. code-block:: bash

        python setup.py clean

        python setup.py install

- run the test  (optional, **run all the examples outside of the installation dir**)

  .. code-block:: bash

       cd ~/
       mkdir test_jetset
       cd test_jetset
       pytest  --pyargs  -vvv jetset.tests.test_users::TestUser





To install from source a C compiler is also necessary, plus the SWIG wrapper generator.

All the dependencies are installed following the Anaconda method **OR** the pip method, as described below.
