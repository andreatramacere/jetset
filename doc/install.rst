.. _install_file:

Installation
============

.. important::
    Starting from version 1.1.0, python 2 is not supported anymore. Python >=3.9 is suggested.


Windows prerequisites
------------------------
Install the Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl/install


Install JetSeT from Anaconda
-----------------------------
to get anaconda: https://www.anaconda.com/download/

or you can use mamba: https://mamba.readthedocs.io/en/latest/

if you use mamba replace ``conda`` with ``mamba`` in the following

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

      conda create --name jetset python=3.10 ipython jupyter

  .. code-block:: bash

      conda activate jetset

- install the code

  .. code-block:: bash

      conda install -c andreatramacere -c astropy -c conda-forge 'jetset>=1.3'



- run the test (optional)

  .. code-block:: bash

    pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser


Install  JetSeT from pip
------------------------------------------------------------------------------

- create a virtual environment (not necessary, but suggested):

  .. code-block:: bash

    pip install virtualenv

    virtualenv -p python3.10 jetset

    source jetset/bin/activate

    pip install ipython jupyter

- MacOS
  
  .. code-block:: bash

      pip install jetset

- Linux
  
  .. code-block:: bash
    
      pip install jetset

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

Download the code
^^^^^^^^^^^^^^^^^

To install from source a C compiler is also necessary, plus the SWIG wrapper generator.

- Get the source code from: https://github.com/andreatramacere/jetset/archive/stable.tar.gz

- Uncompress the  archive:  `jetset-stable.tar.gz`

- cd to  the dir source code dir

  .. code-block:: bash

      cd jetset-stable

Installation from source using Anaconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line:


  .. code-block:: bash

      conda install -c astropy -c conda-forge --file requirements.txt

.. important::
    if anaconda fails to install swig, you can try one of the following alternative :ref:`swig` 


-  run on the command line

   .. code-block:: bash

       pip install .

- run the test (optional, **run all the examples outside the installation directory**)

  .. code-block:: bash

     cd ~/

     mkdir test_jetset

     cd test_jetset

     pytest --disable-warnings  --pyargs  -vvv jetset.tests.test_users::TestUser






Installation from source using PIP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Install requirements, run on the command line: 

  .. code-block:: bash

    pip install -r requirements.txt

.. important::
    if pip fails to install swig, you can try one of the following alternative :ref:`swig` 


- Install JetSeT: run on the command line:

  .. code-block:: bash

        pip install .

- run the test  (optional, **run all the examples outside the installation directory **)

  .. code-block:: bash

       cd ~/
       mkdir test_jetset
       cd test_jetset
       pytest  --pyargs  -vvv jetset.tests.test_users::TestUser







