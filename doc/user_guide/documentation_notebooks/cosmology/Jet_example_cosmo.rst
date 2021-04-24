.. _cosmology:

Choosing a cosmology model
==========================

Since the new implementation allows to switch on the fly from the
``src`` to the ``obs`` frame, we need to define a cosmolgoy model in
order to evaluate the luminosity distance also for dataset independently
from the Jet model. This is implemented in usign the astropy cosmology
subpackaged https://docs.astropy.org/en/stable/cosmology/ Moreover, this
allows also to change easilyt the cosmology used for the Jet model,
passing a specifc cosmology model form astropy

.. code:: ipython3

    import jetset
    print('tested on jetset',jetset.__version__)


.. parsed-literal::

    tested on jetset 1.2.0rc4


.. code:: ipython3

    from jetset.cosmo_tools import Cosmo
    c=Cosmo()

the defualt comoslogy is

.. code:: ipython3

    c




.. parsed-literal::

    FlatLambdaCDM(name="Planck13", H0=67.8 km / (Mpc s), Om0=0.307, Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV, Ob0=0.0483)



you can choose a different cosmology model from astropy eg

.. code:: ipython3

    from astropy.cosmology import WMAP9 as cw
    c=Cosmo(astropy_cosmo=cw)

.. code:: ipython3

    c




.. parsed-literal::

    FlatLambdaCDM(name="WMAP9", H0=69.3 km / (Mpc s), Om0=0.286, Tcmb0=2.725 K, Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.0463)



or you can set directly a luminosity distance in cm

.. code:: ipython3

    c=Cosmo(DL_cm=1e28)


.. parsed-literal::

    using cosmo without z and only DL, should be used only for galactic objects!!
    z will be fixed to zero


.. code:: ipython3

    c




.. parsed-literal::

    cosmology is not defined, the luminosity distance has been set to 1e+28 cm



IF you don’t choose a cosmology model, the default one will be used
(``FlatLambdaCDM(name="Planck13", H0=67.8 km / (Mpc s), Om0=0.307, Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV, Ob0=0.0483)``).

Changing cosmology model in data
--------------------------------

.. code:: ipython3

    from jetset.data_loader import Data
    from jetset.data_loader import ObsData
    from jetset.cosmo_tools import Cosmo
    from jetset.test_data_helper import  test_SEDs
    c=Cosmo()
    data_table=Data.from_file(test_SEDs[1])
    sed_data=ObsData(data_table=data_table,cosmo=c)

Changing cosmology model in Jet models
--------------------------------------

.. code:: ipython3

    from jetset.jet_model import Jet
    my_jet=Jet(cosmo=c)

or for already built models

.. code:: ipython3

    my_jet.cosmo=c
