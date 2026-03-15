.. _install_file:

Installation
============

.. important::
    Python >=3.11,<=3.13 is suggested. Python 3.14 has been tested, no binaries provided, you can install from sources. 


Windows prerequisites
------------------------
Install the Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl/install


Install JetSeT from Anaconda
-----------------------------
.. note:: anaconda binaries for apple M processors are distributed only for python version >= 3.11

to get anaconda: https://www.anaconda.com/download/

or you can use mamba: https://mamba.readthedocs.io/en/latest/

if you use mamba replace ``conda`` with ``mamba`` in the following

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

      conda create --name jetset python=3.12 ipython jupyter

  .. code-block:: bash

      conda activate jetset

- install the code

  .. code-block:: bash

      conda install -c andreatramacere -c astropy -c conda-forge 'jetset>=1.4'



- run the test (optional)

  .. code-block:: bash

    pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser


Install  JetSeT from pip
------------------------------------------------------------------------------

.. note:: pip binaries for apple M processors are distributed only for python version >= 3.11

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

    pip install virtualenv

    virtualenv -p python3.10 jetset

    source jetset/bin/activate

    pip install ipython jupyter

- install the code
  
  .. code-block:: bash

 

  if fails, use one of the following methods 

  - Use anaconda

  OR

  - Install from source
  
  

- run the test (optional)

  .. code-block:: bash

    pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser


Install binaries from GitHub
------------------------------
To use the git release binaries, follow the instructions here: :ref:`install_pre_file`


Install the JetSeT from source
------------------------------
As for the other cases, I suggest creating a specific virtual environment

Download the code
^^^^^^^^^^^^^^^^^

To install from source a C compiler is also necessary, plus the SWIG wrapper generator (the latter will be installed by dependencies).

- Typically, you will clone the repo (and checkout a specific branch if needed).  In case you want access a specific release/tag e.g. : https://github.com/andreatramacere/jetset/archive/refs/tags/stable.tar.gz

- Uncompress the  archive:  `jetset-stable.tar.gz`

- cd to  source code dir

  .. code-block:: bash

      cd jetset-stable

Installation from source using the `install.sh` script (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- run on the command line:


  .. code-block:: bash

      ./install.sh

The `install.sh` script will take care of everything from installing dependencies to compile and install jetset.


- run the test (optional, **run all the examples outside the installation directory**)

  .. code-block:: bash

     cd ~/

     mkdir test_jetset

     cd test_jetset

     pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser

.. important::
    if pip or conda fails to install swig, you can try one of the following alternative :ref:`swig` 


Installation from source without using the `install.sh` script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- install the dependencies:

  - for pip: ``pip install -r requirements.txt``

  OR

  - for conda: ``conda install -c astropy -c conda-forge --file requirements.txt``


- install the code:

  ``pip install .``




