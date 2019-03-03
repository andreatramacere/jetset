
.. _data-format:

Data format and SED data
========================

The data are treated using two classes from the module :mod:`.data_loader`

- :class:`jetset.data_loader.Data` class 
- :class:`jetset.data_loader.ObsData` class 

The class :class:`jetset.data_loader.Data` is in charge of storing the data, giving access to the I/O fucntionalities, and provides and interface  to the `~astropy.table.Table` (see the `astropy`_ documentation, for further information)

The class :class:`jetset.data_loader.ObsData` use the information stored in :class:`jetset.data_loader.Data`, and can perform several operations 
 
 - rebinning of the data
 - selection of time ranges
 - selection of datasets
 - transformation from liner to logarthmic representation
 - handling of errors and systematcis
 

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')
    
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline  

data format
-----------

The SED data are internally stored as astropy tables, but it is very
easy to import from

1. ascii files
2. numpy array in general

once that is clear the data format. The easiest way to undertand the
data format is to build an empty table to have a look at the structure
of the table:

.. code:: ipython3

    from jetset.data_loader import Data
    data=Data(n_rows=10)

we can easily access the astropy table

.. code:: ipython3

    data.table




.. raw:: html

    <i>Table length=10</i>
    <table id="table4610276880" class="table-striped table-bordered table-condensed">
    <thead><tr><th>x</th><th>dx</th><th>y</th><th>dy</th><th>T_start</th><th>T_stop</th><th>UL</th><th>data_set</th></tr></thead>
    <thead><tr><th>Hz</th><th>erg / (cm2 s)</th><th>Hz</th><th>erg / (cm2 s)</th><th>MJD</th><th>MJD</th><th></th><th></th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>bool</th><th>bytes16</th></tr></thead>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    <tr><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>False</td><td>0.0</td></tr>
    </table>



-  ``x`` column is reserved to frequencies (mandatory)
-  ``y`` columm is reserved to fluxes (mandatory)
-  ``dx`` columm is reserved to the error on the frequency,or bin width
-  ``dy`` columm is reserved to the error on the fluxes
-  ``UL`` columm is reserved to the flag for Upper Limit
-  ``T_start`` and ``T_stop`` are used to indentify the time range to
   select data using the class ``ObsData``
-  ``data_set``

and we can easily access the metadata

.. code:: ipython3

    data.metadata




.. parsed-literal::

    OrderedDict([('z', 0),
                 ('UL_CL', 0.95),
                 ('restframe', 'obs'),
                 ('data_scale', 'lin-lin'),
                 ('obj_name', 'new-src')])



-  ``z``: the redshift of the object
-  ``UL_CL``: the CL for the UL
-  ``restframe``: possible values\ ``obs`` or ``src``, indicating if the
   data are oserved flux, or luminosities, respectively
-  ``data_scale``: ossible values\ ``lin-lin`` or ``log-log``,
   indicating if the data are in linear or logarithmic scale,
   respectively
-  ``obj_name``: the name of the object

**we remind that the conversion from ``src`` to ``obs`` will be exposed
in the next release, for the time being, please use only fluxes and not
luminosities**

loading from astropy table
--------------------------

.. code:: ipython3

    from jetset.test_data_helper import  test_SEDs
    test_SEDs




.. parsed-literal::

    ['/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-1.0.2-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_3C345.dat',
     '/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-1.0.2-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk421.dat',
     '/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-1.0.2-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk501.dat']



As you can see there are three 3 files. We use in our example the file for Mrk 421, and we use class:`jetset.data_loader.Data` class  

.. code:: ipython3

    from jetset.data_loader import Data

.. code:: ipython3

    data=Data(data_table=test_SEDs[1])

.. code:: ipython3

    data.table




.. raw:: html

    <i>Table length=110</i>
    <table id="table90555021520" class="table-striped table-bordered table-condensed">
    <thead><tr><th>x</th><th>dx</th><th>y</th><th>dy</th><th>T_start</th><th>T_stop</th><th>UL</th><th>data_set</th></tr></thead>
    <thead><tr><th>Hz</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>MJD</th><th>MJD</th><th></th><th></th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>bool</th><th>str13</th></tr></thead>
    <tr><td>2299540000.0</td><td>0.0</td><td>1.3409e-14</td><td>3.91e-16</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>2639697000.0</td><td>0.0</td><td>1.793088e-14</td><td>3.231099e-26</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>4799040000.0</td><td>0.0</td><td>2.3136e-14</td><td>2.4e-16</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>4805039000.0</td><td>0.0</td><td>1.773414e-14</td><td>1.773414e-15</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>4843552000.0</td><td>0.0</td><td>2.77614e-14</td><td>2.615339e-26</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>7698460000.0</td><td>0.0</td><td>3.696e-14</td><td>4.62e-16</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>8267346000.0</td><td>0.0</td><td>2.836267e-14</td><td>2.836267e-15</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>8331867000.0</td><td>0.0</td><td>3.98963e-14</td><td>3.627671e-26</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>8388659000.0</td><td>0.0</td><td>3.16345e-14</td><td>1.931495e-15</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>2.417992e+25</td><td>0.0</td><td>9.754259e-11</td><td>3.560456e-11</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>3.823193e+25</td><td>0.0</td><td>8.199207e-11</td><td>7.050657e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>6.059363e+25</td><td>0.0</td><td>5.614334e-11</td><td>5.793969e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>6.073707e+25</td><td>0.0</td><td>1.14705e-10</td><td>6.573696e-11</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>9.603433e+25</td><td>0.0</td><td>4.662219e-11</td><td>5.097912e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>1.522041e+26</td><td>0.0</td><td>5.221583e-11</td><td>4.89063e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>2.41227e+26</td><td>0.0</td><td>3.66834e-11</td><td>4.682033e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>3.823193e+26</td><td>0.0</td><td>2.247871e-11</td><td>4.343216e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>6.059363e+26</td><td>0.0</td><td>1.972081e-11</td><td>4.407365e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    <tr><td>9.603433e+26</td><td>0.0</td><td>7.994215e-12</td><td>3.469109e-12</td><td>0.0</td><td>0.0</td><td>False</td><td>campaing-2009</td></tr>
    </table>



.. code:: ipython3

    data.metadata




.. parsed-literal::

    OrderedDict([('z', 0.0308),
                 ('restframe', 'obs'),
                 ('data_scale', 'lin-lin'),
                 ('obj_name', 'J1104+3812,Mrk421')])



this is an extract of the  astropy table saved in the format ``ascii.ecsv`` 

::

    # %ECSV 0.9
    # ---
    # datatype:
    # - {name: x, unit: Hz, datatype: float64}
    # - {name: y, unit: erg/(cm2 s) , datatype: float64}
    # - {name: dy, unit: erg/(cm2 s) ,datatype: float64}
    # - {name: data_set, datatype: string}
    # meta: !!omap
    # - {z: 0.0308 }
    # - {restframe: obs}
    # - {data_scale: lin-lin}
    # - {obj_name: 'J1104+3812,Mrk421'}
    # schema: astropy-2.0
    x y dy data_set
    2.299540e+09 1.340900e-14 3.910000e-16    campaing-2009
    2.639697e+09 1.793088e-14 3.231099e-26    campaing-2009
    4.799040e+09 2.313600e-14 2.400000e-16    campaing-2009

saving to a file
----------------

importing data from an arbitrary ascii file or numpy array
----------------------------------------------------------

Assume that your data are stored in an ASCII file named
'test-ascii.txt', with - ``x`` in the first colum with frequency in
``Hz`` , - ``y`` in the second column with fluxes in erf ``cm-2 s-1``, -
``dy`` in the third column with the same units as ``y`` - the data are
in ``log-log`` scale

**of course the column number depends on the file that you are using,
this is only an example**

.. code:: ipython3

    from jetset.data_loader import Data
    import numpy as np
    
    d=np.genfromtxt('test-ascii.txt')
    data=Data(n_rows=d.shape[0])
    data.set_field('x',d[:,0])
    data.set_field('y',d[:,1])
    data.set_field('dy',value=d[:,2])

then you can set the meatdata as follows

.. code:: ipython3

    data.set_meta_data('z',1.02)
    data.set_meta_data('restframe','obs')
    data.set_meta_data('data_scale','log-log')


of course this method applies if you have a generic 2-dim numpy array.

importing data from the asi ssdc sedtool
----------------------------------------

To import data from a data file downloaded from the asi ssdc sedtool:
https://tools.ssdc.asi.it/SED/

we can use the importing tool in the :class:`jetset.data_loader.Data`. We just need to have the file downloaded from the asi ssdc sedtool, and to know the redshfit of the object, the scale we selected (lin-lin, or log-log).
Assume that we downloaded the data for Mrk421, in observed fluxes and linear scale, and the data are saved in the file 'MRK421_asdc.txt', we only have to do:

.. code:: ipython3

    from jetset.data_loader import Data
    data=Data.from_asdc(asdc_sed_file='MRK421_asdc.txt',obj_name='Mrk421',restframe='obs',data_scale='lin-lin',z=0.038)


**we remind that the conversion from ``src`` to ``obs`` will be exposed
in the next release, for the time being, please use only fluxes and not
luminosities**

Building the SED
----------------

Once we have a data table built with the class:`jetset.data_loader.Data`, following  one of the method described above, you can create  SED data using the  :class:`jetset.data_loader.ObsData` class.
In the example we use one of the test SEDs provided by the package:

We start to loading the SED of Mrk 421, and we pass to ``ObsData``
directly the path to the file, because this is already in the format
that we need and that we have discussed above.

.. code:: ipython3

    from jetset.data_loader import Data
    from jetset.data_loader import ObsData
    from jetset.test_data_helper import  test_SEDs
    
    data_table=Data(data_table=test_SEDs[1])
    sed_data=ObsData(data_table=data_table)

As you can see the all the meta-data have been properly sourced from the
SED file header. You also get information on the lenght of the data,
before and after elimination of duplicated entries, and upper limits
These meta-data are parameters needed by the

.. code:: ipython3

    sed_data.table




.. raw:: html

    <i>Table length=110</i>
    <table id="table90560932944" class="table-striped table-bordered table-condensed">
    <thead><tr><th>nu_data</th><th>dnu_data</th><th>nuFnu_data</th><th>dnuFnu_data</th><th>nu_data_log</th><th>dnu_data_log</th><th>nuFnu_data_log</th><th>dnuFnu_data_log</th><th>dnuFnu_facke</th><th>dnuFnu_facke_log</th><th>UL</th><th>zero_error</th><th>T_start</th><th>T_stop</th><th>data_set</th></tr></thead>
    <thead><tr><th>Hz</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>Hz</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th></th><th></th><th></th><th>MJD</th><th>MJD</th><th></th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>bool</th><th>bool</th><th>float64</th><th>float64</th><th>str13</th></tr></thead>
    <tr><td>2299540000.0</td><td>0.0</td><td>1.3409e-14</td><td>3.91e-16</td><td>9.361640968434164</td><td>0.0</td><td>-13.872603609223393</td><td>0.012663818511758627</td><td>2.6818000000000003e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>2639697000.0</td><td>0.0</td><td>1.793088e-14</td><td>3.231099e-26</td><td>9.421554078847052</td><td>0.0</td><td>-13.746398395894273</td><td>7.825876176646739e-13</td><td>3.586176e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>4799040000.0</td><td>0.0</td><td>2.3136e-14</td><td>2.4e-16</td><td>9.681154369792159</td><td>0.0</td><td>-13.635711724385564</td><td>0.0045051294803241885</td><td>4.627200000000001e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>4805039000.0</td><td>0.0</td><td>1.773414e-14</td><td>1.773414e-15</td><td>9.68169691696108</td><td>0.0</td><td>-13.751189867373059</td><td>0.04342944819032518</td><td>3.546828e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>4843552000.0</td><td>0.0</td><td>2.77614e-14</td><td>2.615339e-26</td><td>9.68516396664987</td><td>0.0</td><td>-13.556558636309997</td><td>4.091390549490907e-13</td><td>5.55228e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>7698460000.0</td><td>0.0</td><td>3.696e-14</td><td>4.62e-16</td><td>9.886403857589054</td><td>0.0</td><td>-13.43226803745193</td><td>0.005428681023790648</td><td>7.392e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>8267346000.0</td><td>0.0</td><td>2.836267e-14</td><td>2.836267e-15</td><td>9.917366113839973</td><td>0.0</td><td>-13.547252888027566</td><td>0.043429448190325175</td><td>5.672534000000001e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>8331867000.0</td><td>0.0</td><td>3.98963e-14</td><td>3.627671e-26</td><td>9.920742328771254</td><td>0.0</td><td>-13.399067379102538</td><td>3.948931348171262e-13</td><td>7.97926e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>8388659000.0</td><td>0.0</td><td>3.16345e-14</td><td>1.931495e-15</td><td>9.92369254063231</td><td>0.0</td><td>-13.499839025404517</td><td>0.026516544289422034</td><td>6.3268999999999995e-15</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>2.417992e+25</td><td>0.0</td><td>9.754259e-11</td><td>3.560456e-11</td><td>25.38345485965064</td><td>0.0</td><td>-10.010805716985434</td><td>0.15852422965797036</td><td>1.9508518000000003e-11</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>3.823193e+25</td><td>0.0</td><td>8.199207e-11</td><td>7.050657e-12</td><td>25.582426222350527</td><td>0.0</td><td>-10.086228149101405</td><td>0.03734582416192853</td><td>1.6398414000000002e-11</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>6.059363e+25</td><td>0.0</td><td>5.614334e-11</td><td>5.793969e-12</td><td>25.78242697068017</td><td>0.0</td><td>-10.250701754501332</td><td>0.044819007294872405</td><td>1.1228668000000001e-11</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>6.073707e+25</td><td>0.0</td><td>1.14705e-10</td><td>6.573696e-11</td><td>25.78345383740898</td><td>0.0</td><td>-9.94041765075539</td><td>0.24889236724724106</td><td>2.2941000000000003e-11</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>9.603433e+25</td><td>0.0</td><td>4.662219e-11</td><td>5.097912e-12</td><td>25.982426510793527</td><td>0.0</td><td>-10.33140733007377</td><td>0.04748801055523926</td><td>9.324438000000001e-12</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>1.522041e+26</td><td>0.0</td><td>5.221583e-11</td><td>4.89063e-12</td><td>26.1824263514056</td><td>0.0</td><td>-10.282197814249994</td><td>0.04067681433064456</td><td>1.0443166e-11</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>2.41227e+26</td><td>0.0</td><td>3.66834e-11</td><td>4.682033e-12</td><td>26.38242591580127</td><td>0.0</td><td>-10.43553041856344</td><td>0.05543055158433863</td><td>7.33668e-12</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>3.823193e+26</td><td>0.0</td><td>2.247871e-11</td><td>4.343216e-12</td><td>26.582426222350527</td><td>0.0</td><td>-10.648228615520983</td><td>0.08391205467368516</td><td>4.495742000000001e-12</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>6.059363e+26</td><td>0.0</td><td>1.972081e-11</td><td>4.407365e-12</td><td>26.78242697068017</td><td>0.0</td><td>-10.705075251093293</td><td>0.09705961870904517</td><td>3.944162e-12</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    <tr><td>9.603433e+26</td><td>0.0</td><td>7.994215e-12</td><td>3.469109e-12</td><td>26.982426510793527</td><td>0.0</td><td>-11.097224175808465</td><td>0.18846314438889974</td><td>1.598843e-12</td><td>0.2</td><td>False</td><td>False</td><td>0.0</td><td>0.0</td><td>campaing-2009</td></tr>
    </table>



.. code:: ipython3

    sed_data.metadata


.. parsed-literal::

    z :  0.0308
    obj_name :  J1104+3812,Mrk421
    restframe :  obs
    data_scale :  lin-lin
    UL_CL :  0.95


Plotting data
-------------

We can now plot our SED using the :class:`BlazarSEDFit.plot_sedfit.Plot` class 


.. code:: ipython3

    from jetset.plot_sedfit import PlotSED
    myPlot=PlotSED(sed_data)



.. image:: Jet_example_load_data_files/Jet_example_load_data_43_0.png


or you can create the object to plot on the fly in this way

.. code:: ipython3

    myPlot=sed_data.plot_sed()




.. image:: Jet_example_load_data_files/Jet_example_load_data_45_0.png


you can rescale your plot

.. code:: ipython3

    myPlot=sed_data.plot_sed()
    myPlot.rescale(x_min=6,x_max=28,y_min=-16,y_max=-9)



.. image:: Jet_example_load_data_files/Jet_example_load_data_47_0.png


**to have interactive plot in jupyter**

.. code:: ipython3

    %matplotlib notebook
    myPlot=sed_data.plot_sed()



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAtAAAAGwCAYAAACAS1JbAAAgAElEQVR4Xu2dB7hdRbn+35MeCBACwUAkBJUk1CQE6QIBFZAiSlOke6XIxYL04o3XiwGBK4jAJUpT/gi5gFIE6UWkXAmEhBaQJgIiLUBMSP8/X1wHdzb77FVmrTWz1vnt5zkP4eyZ7/vW75195j1zZs3qEC8IQAACEIAABCAAAQhAIDGBjsQtaQgBCEAAAhCAAAQgAAEICAPNIIAABCAAAQhAAAIQgEAKAhjoFLBoCgEIQAACEIAABCAAAQw0YwACEIAABCAAAQhAAAIpCGCgU8CiKQQgAAEIQAACEIAABDDQjAEIQAACEIAABCAAAQikIICBTgGLphCAAAQgAAEIQAACEMBAMwYgAAEIQAACEIAABCCQggAGOgUsmkIAAhCAAAQgAAEIQAADzRiAAAQgAAEIQAACEIBACgIY6BSwaAoBCEAAAhCAAAQgAAEMNGMAAhCAAAQgAAEIQAACKQhgoFPAoikEIAABCEAAAhCAAAQw0IwBCEAAAhCAAAQgAAEIpCCAgU4Bi6YQgAAEIAABCEAAAhDAQDMGIAABCEAAAhCAAAQgkIIABjoFLJpCAAIQgAAEIAABCEAAA80YgAAEIAABCEAAAhCAQAoCGOgUsGgKAQhAAAIQgAAEIAABDDRjAAIQgAAEIAABCEAAAikIYKBTwKIpBCAAAQhAAAIQgAAEMNCMAQhAoHQCTzzxRJ958+b9vKOjY0tJPUsvgIQQWJrAwsWLF9/Xp0+fb6y77rrzgAMBCEAgjgAGOo4Q70MAArkTePTRR789cODAI4cNG/Zujx49FueegIAQSEFg0aJFHS+99NLAd99996djx449J0VXmkIAAt2UAAa6mwrPZUPAJ4GpU6dOHTVqVEffvn3n+6yD3BDoJDB37tzeM2bMWDR69OixUIEABCAQRwADHUeI9yEAgdwJTJ069YUNNtjgrY4OfgTlDpeAmQgsXrxY06ZNGzRmzJhPZApAJwhAoFsRYPbqVnJzsRAIg8DUqVNfHD169JthVEMVEPgngccee2zlMWPGDIcHBCAAgTgCGOg4QrwPAQjkTsDVQH/xZ/eNtKKu+/ctZ+ReHAFLJ/DTn/50pcsuu2zlKVOmeNUTA1269CSEQGUJYKArKx2FQ6C6BEI30JtsssnIvffe+60jjjjirS9+8YtrTp8+fdlXX321zw033PDMzjvv/H4n+euvv365U089ddUnn3xymeWXX37hK6+8Mr1RlRkzZvTZf//9h0+bNm3ZIUOGzPvJT37yl9122+3D/p1tN9100xEPPfTQcvPmzZvSu3fvVMJusskmI5555pn+8+fP7zF06NC5p5xyyqv77rvvzM4gp5566irnn3/+x959991ew4cP/+Css856efvtt59l77er/5VXXul12GGHrW51zZkzp8daa631wZlnnvnytttu+49WBe6+++7Dr7322pUuv/zyP3/ta197t7PNwQcfvPoll1yyyjnnnPPit771rbda9U1qoPfYY4/h11xzzUrTp09/fL311ps7Z86cjgMOOGDYfffdt7xd37Bhw+ZOmDDhr3vvvfd7lmfu3Lkd7fRrrgUDnWro0RgC3ZoABrpby8/FQ8APgSoZ6NNPP33wJptsMnvffff9xCWXXPJCo4G+8847l3nqqaf6mcH8yU9+smqzgR4zZsyojTbaaNbZZ5/9ytVXX73CkUceOXzGjBmPr7baags6yV9wwQWDLrroosFTpkwZ0JWBNnO69dZbv9/KgD700EP9N9xwwzlmvO+8885ld9lllxFPPvnk42usscZ8+/+ddtppxK233jpjiy22mH3GGWcMPu2001Z74403HuvVq5e177L+J598ss/kyZNXPOigg94eOnTo/LPPPnvlU089deiLL744fYUVVljUPHKsRruGkSNHzrnlllues/fnz5+voUOHbtCvX79FRx999Gut6rc2F1xwQewK9C233DLgpJNOGmo5Og30e++912PChAlDDj300DfXWmuteZMnT17h61//+iceeeSRJ0aOHDnPDHQ7/TDQfj7/ZIVAHQhgoOugItcAgYoRcDHQlz/40qAf3PDE8PkLF3esslzfed/abq1X9t10jbfzRNC5An3UUUd9uE/7Yx/72AYXXXTRUga6M+dvf/vb5Y444ojhjQZ62rRpfTfaaKN1X3/99akrrrjiEsM5btw4W9l++9hjj33D/v/tt9/uueGGG6598cUXv7DddtuNymKgG6/bDPGOO+446uabb3562223nT1p0qQVzz333CHTp09/ytqZ4VxhhRXGvvjii9PMYLervxXPAQMGjL355ptnfOYzn5ndykCvtNJKC2wV+qmnnnp88ODBC6+88soVzj///FX+8Y9/9DjggAPeNANtq82XXnrp4LFjx/7j6quvXumAAw74+6c+9am5jVs4Dj300I8//PDDy952221/HjRo0EIz2RtssME6l1566QubbrrpOp0GulWNI0aMWOfEE0989cADD/xwFd7atdOvMw4r0Hl+iogFgXoTwEDXW1+uDgJBEshqoM08//DGJ9eYu2BRj84L69urx6JTdl7npTxNdB4G+pe//OXACRMmDH3++eef6Kx1//33H9bR0bH4sssue9m+t99++w375Cc/+cFXv/rVmaNGjVo/q4EeP378p+6///7l582b17Hlllu+d8899zzbo0cPvfXWWz222mqrkeeee+5fttpqq3+cfvrpq1x++eUrP/HEE0/a+2kM9AMPPNB/m222Wfu11157zExtKwM9dOjQeW+88UbvMWPGzD7uuOPe+MIXvvCJL37xizMnTZo0uNFAH3XUUcN/8IMfvHzsscf+3VaJL7744kFmoP/0pz/N2Geffdb461//2uemm256bvnll1/yi8cpp5zyMfveJZdc8nJHR8e4rgz0yy+/3OtTn/rUBg8++OCTY8eO/aCxRgx0kD8KKAoClSWAga6sdBQOgeoSSGOgO28YtKt98rX3lrWV5+Yr792zY/E6qy6/ZG9uHjcW5mGgzzvvvEGTJk1a5bHHHnu6s94jjzxy6Kuvvtr7mmuuefHee+9d5pBDDhk+ffr0J59//vk+Lgba4psRtZVw21IyYcKEv9v3Fi1apBNPPHHImWeeuZqkjuWWW27Bb3/722e33nrrpVaQW62gNzI2I77ZZpuN2n333d+eOHHi31qNPNvCYQZ6p512eu+44477+J133vnsWmuttd5LL700beuttx7RaKAnTpy42muvvfbhfnFblf7FL34xePXVV5+3YMGCjuuvv/75vn37LnnAzrPPPtt7u+22Gzl16tSnzLh3ZaDt+sePH7/W8OHD515xxRUvNdeIga7uzwsqh0CIBDDQIapCTRCoOYGsBvqxv747oCs0oz++wpIb40Ix0LYC/YMf/GDoc8899+EK9AEHHLC61WgrqaNHj1574sSJL++8886z7GbDZgNtWxFee+21Ptb+gw8+6NGrV6/F9mX/v9tuu739q1/96i+tWHzmM59Z67DDDvu73ch31llnrXzOOecMufHGG5+1m+6uvfba5Q899NA1p0yZ8uTw4cMTbeGYNWuWGdMRtlJ+5ZVXfsSYdtbQaaB/+tOfvjps2LD1dtxxx5nvvPNOz8mTJ79kW1eat3A88sgjH/5iYQb6+9///uq2l/zee+99apNNNpnTGXf77bf/5K677jrTbui077Uy0PaLwq677vqJ999/v8ett976XKf5buSDga75DxUuDwIlE8BAlwycdBCAgJTGQDfy2vjU29f/+/tzl5jKxpfthf6/kz671AkYLpzzWIG2PdCf/vSn1/3b3/724R7ojTbaaORee+319sEHH/z2KqusMmbQoEFLbiZcuHChZs6c2cv2EF9++eXP7bDDDkt+GWg0p13dRNh8nZtvvvmIz3/+8zNtFdq2jPTu3XvxRRddtGTLiL1GjRq1znHHHffaQQcd9E7n97pagbZTLj73uc99asUVV1xw3XXXvdC47aM5b6OBPuqoo1Y7++yzV73++utn2C8IzQa6+ci6zlM4DjnkkDdOOeWUj992220zRo8ePddyLLfccmP69OmzuPOhO2+99VavgQMHLvjRj3708uGHH/62mee99tpr+Msvv9znjjvueHbAgAEtHw2PgXb5RNAXAhBoJoCBZkxAAAKlE8hqoH3sgTYTaU+pW3PNNde/4IILXtxhhx3e79ev32Izk2bePvjgg44bbrhhue9+97tr/PnPf368Z8+e6lwBHT169KhNNtnkw1M47EZDO4VjyJAhC+yYuE7wtoXD9hc///zz0+yEjuYV1K5O4Xj00Uf7PfPMM3122mmn9+0UjosuumjFb3/728Pvuuuup7fccsvZ55577kpnnXXWqjfddNMzo0aNmnfdddctv88++3zy/vvvf8r2CLer37ZE7Ljjjp/s0aPH4ptvvvm55uP1OlfNn3766el24kWjgX799dd7Pvjgg8vssssu7xunpAbazoG2mu20j9tvv32GrZobJ6uz8zVs2LDRd9xxx9Mbb7zxbDPL++yzz7AnnnhimXvvvfeZVqeDtNOveeBzE2HpPwpICIHKEsBAV1Y6CodAdQlkNdB2xWWdwvGVr3zlre9+97tvDh06dH07A7qRdqdpvPHGG5ezY+Ma39t4441nPfTQQ0seCBKdA73mY489tuQc6LPPPvulVudAt9rC0RizKwM9ZcqUfgcddNDw5557rr8Z3TXWWGPuscce+9r++++/5AQKM562GnzVVVet9N57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" width="720">


**to have interactive plot in ipython**

.. code:: ipython3

    from matplotlib import pylab as plt
    plt.ion()
    myPlot=sed_data.plot_sed()



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAtAAAAGwCAYAAACAS1JbAAAgAElEQVR4Xu2dB7hdRbn+35MeCBACwUAkBJUk1CQE6QIBFZAiSlOke6XIxYL04o3XiwGBK4jAJUpT/gi5gFIE6UWkXAmEhBaQJgIiLUBMSP8/X1wHdzb77FVmrTWz1vnt5zkP4eyZ7/vW75195j1zZs3qEC8IQAACEIAABCAAAQhAIDGBjsQtaQgBCEAAAhCAAAQgAAEICAPNIIAABCAAAQhAAAIQgEAKAhjoFLBoCgEIQAACEIAABCAAAQw0YwACEIAABCAAAQhAAAIpCGCgU8CiKQQgAAEIQAACEIAABDDQjAEIQAACEIAABCAAAQikIICBTgGLphCAAAQgAAEIQAACEMBAMwYgAAEIQAACEIAABCCQggAGOgUsmkIAAhCAAAQgAAEIQAADzRiAAAQgAAEIQAACEIBACgIY6BSwaAoBCEAAAhCAAAQgAAEMNGMAAhCAAAQgAAEIQAACKQhgoFPAoikEIAABCEAAAhCAAAQw0IwBCEAAAhCAAAQgAAEIpCCAgU4Bi6YQgAAEIAABCEAAAhDAQDMGIAABCEAAAhCAAAQgkIIABjoFLJpCAAIQgAAEIAABCEAAA80YgAAEIAABCEAAAhCAQAoCGOgUsGgKAQhAAAIQgAAEIAABDDRjAAIQgAAEIAABCEAAAikIYKBTwKIpBCAAAQhAAAIQgAAEMNCMAQhAoHQCTzzxRJ958+b9vKOjY0tJPUsvgIQQWJrAwsWLF9/Xp0+fb6y77rrzgAMBCEAgjgAGOo4Q70MAArkTePTRR789cODAI4cNG/Zujx49FueegIAQSEFg0aJFHS+99NLAd99996djx449J0VXmkIAAt2UAAa6mwrPZUPAJ4GpU6dOHTVqVEffvn3n+6yD3BDoJDB37tzeM2bMWDR69OixUIEABCAQRwADHUeI9yEAgdwJTJ069YUNNtjgrY4OfgTlDpeAmQgsXrxY06ZNGzRmzJhPZApAJwhAoFsRYPbqVnJzsRAIg8DUqVNfHD169JthVEMVEPgngccee2zlMWPGDIcHBCAAgTgCGOg4QrwPAQjkTsDVQH/xZ/eNtKKu+/ctZ+ReHAFLJ/DTn/50pcsuu2zlKVOmeNUTA1269CSEQGUJYKArKx2FQ6C6BEI30JtsssnIvffe+60jjjjirS9+8YtrTp8+fdlXX321zw033PDMzjvv/H4n+euvv365U089ddUnn3xymeWXX37hK6+8Mr1RlRkzZvTZf//9h0+bNm3ZIUOGzPvJT37yl9122+3D/p1tN9100xEPPfTQcvPmzZvSu3fvVMJusskmI5555pn+8+fP7zF06NC5p5xyyqv77rvvzM4gp5566irnn3/+x959991ew4cP/+Css856efvtt59l77er/5VXXul12GGHrW51zZkzp8daa631wZlnnvnytttu+49WBe6+++7Dr7322pUuv/zyP3/ta197t7PNwQcfvPoll1yyyjnnnPPit771rbda9U1qoPfYY4/h11xzzUrTp09/fL311ps7Z86cjgMOOGDYfffdt7xd37Bhw+ZOmDDhr3vvvfd7lmfu3Lkd7fRrrgUDnWro0RgC3ZoABrpby8/FQ8APgSoZ6NNPP33wJptsMnvffff9xCWXXPJCo4G+8847l3nqqaf6mcH8yU9+smqzgR4zZsyojTbaaNbZZ5/9ytVXX73CkUceOXzGjBmPr7baags6yV9wwQWDLrroosFTpkwZ0JWBNnO69dZbv9/KgD700EP9N9xwwzlmvO+8885ld9lllxFPPvnk42usscZ8+/+ddtppxK233jpjiy22mH3GGWcMPu2001Z74403HuvVq5e177L+J598ss/kyZNXPOigg94eOnTo/LPPPnvlU089deiLL744fYUVVljUPHKsRruGkSNHzrnlllues/fnz5+voUOHbtCvX79FRx999Gut6rc2F1xwQewK9C233DLgpJNOGmo5Og30e++912PChAlDDj300DfXWmuteZMnT17h61//+iceeeSRJ0aOHDnPDHQ7/TDQfj7/ZIVAHQhgoOugItcAgYoRcDHQlz/40qAf3PDE8PkLF3esslzfed/abq1X9t10jbfzRNC5An3UUUd9uE/7Yx/72AYXXXTRUga6M+dvf/vb5Y444ojhjQZ62rRpfTfaaKN1X3/99akrrrjiEsM5btw4W9l++9hjj33D/v/tt9/uueGGG6598cUXv7DddtuNymKgG6/bDPGOO+446uabb3562223nT1p0qQVzz333CHTp09/ytqZ4VxhhRXGvvjii9PMYLervxXPAQMGjL355ptnfOYzn5ndykCvtNJKC2wV+qmnnnp88ODBC6+88soVzj///FX+8Y9/9DjggAPeNANtq82XXnrp4LFjx/7j6quvXumAAw74+6c+9am5jVs4Dj300I8//PDDy952221/HjRo0EIz2RtssME6l1566QubbrrpOp0GulWNI0aMWOfEE0989cADD/xwFd7atdOvMw4r0Hl+iogFgXoTwEDXW1+uDgJBEshqoM08//DGJ9eYu2BRj84L69urx6JTdl7npTxNdB4G+pe//OXACRMmDH3++eef6Kx1//33H9bR0bH4sssue9m+t99++w375Cc/+cFXv/rVmaNGjVo/q4EeP378p+6///7l582b17Hlllu+d8899zzbo0cPvfXWWz222mqrkeeee+5fttpqq3+cfvrpq1x++eUrP/HEE0/a+2kM9AMPPNB/m222Wfu11157zExtKwM9dOjQeW+88UbvMWPGzD7uuOPe+MIXvvCJL37xizMnTZo0uNFAH3XUUcN/8IMfvHzsscf+3VaJL7744kFmoP/0pz/N2Geffdb461//2uemm256bvnll1/yi8cpp5zyMfveJZdc8nJHR8e4rgz0yy+/3OtTn/rUBg8++OCTY8eO/aCxRgx0kD8KKAoClSWAga6sdBQOgeoSSGOgO28YtKt98rX3lrWV5+Yr792zY/E6qy6/ZG9uHjcW5mGgzzvvvEGTJk1a5bHHHnu6s94jjzxy6Kuvvtr7mmuuefHee+9d5pBDDhk+ffr0J59//vk+Lgba4psRtZVw21IyYcKEv9v3Fi1apBNPPHHImWeeuZqkjuWWW27Bb3/722e33nrrpVaQW62gNzI2I77ZZpuN2n333d+eOHHi31qNPNvCYQZ6p512eu+44477+J133vnsWmuttd5LL700beuttx7RaKAnTpy42muvvfbhfnFblf7FL34xePXVV5+3YMGCjuuvv/75vn37LnnAzrPPPtt7u+22Gzl16tSnzLh3ZaDt+sePH7/W8OHD515xxRUvNdeIga7uzwsqh0CIBDDQIapCTRCoOYGsBvqxv747oCs0oz++wpIb40Ix0LYC/YMf/GDoc8899+EK9AEHHLC61WgrqaNHj1574sSJL++8886z7GbDZgNtWxFee+21Ptb+gw8+6NGrV6/F9mX/v9tuu739q1/96i+tWHzmM59Z67DDDvu73ch31llnrXzOOecMufHGG5+1m+6uvfba5Q899NA1p0yZ8uTw4cMTbeGYNWuWGdMRtlJ+5ZVXfsSYdtbQaaB/+tOfvjps2LD1dtxxx5nvvPNOz8mTJ79kW1eat3A88sgjH/5iYQb6+9///uq2l/zee+99apNNNpnTGXf77bf/5K677jrTbui077Uy0PaLwq677vqJ999/v8ett976XKf5buSDga75DxUuDwIlE8BAlwycdBCAgJTGQDfy2vjU29f/+/tzl5jKxpfthf6/kz671AkYLpzzWIG2PdCf/vSn1/3b3/724R7ojTbaaORee+319sEHH/z2KqusMmbQoEFLbiZcuHChZs6c2cv2EF9++eXP7bDDDkt+GWg0p13dRNh8nZtvvvmIz3/+8zNtFdq2jPTu3XvxRRddtGTLiL1GjRq1znHHHffaQQcd9E7n97pagbZTLj73uc99asUVV1xw3XXXvdC47aM5b6OBPuqoo1Y7++yzV73++utn2C8IzQa6+ci6zlM4DjnkkDdOOeWUj992220zRo8ePddyLLfccmP69OmzuPOhO2+99VavgQMHLvjRj3708uGHH/62mee99tpr+Msvv9znjjvueHbAgAEtHw2PgXb5RNAXAhBoJoCBZkxAAAKlE8hqoH3sgTYTaU+pW3PNNde/4IILXtxhhx3e79ev32Izk2bePvjgg44bbrhhue9+97tr/PnPf368Z8+e6lwBHT169KhNNtnkw1M47EZDO4VjyJAhC+yYuE7wtoXD9hc///zz0+yEjuYV1K5O4Xj00Uf7PfPMM3122mmn9+0UjosuumjFb3/728Pvuuuup7fccsvZ55577kpnnXXWqjfddNMzo0aNmnfdddctv88++3zy/vvvf8r2CLer37ZE7Ljjjp/s0aPH4ptvvvm55uP1OlfNn3766el24kWjgX799dd7Pvjgg8vssssu7xunpAbazoG2mu20j9tvv32GrZobJ6uz8zVs2LDRd9xxx9Mbb7zxbDPL++yzz7AnnnhimXvvvfeZVqeDtNOveeBzE2HpPwpICIHKEsBAV1Y6CodAdQlkNdB2xWWdwvGVr3zlre9+97tvDh06dH07A7qRdqdpvPHGG5ezY+Ma39t4441nPfTQQ0seCBKdA73mY489tuQc6LPPPvulVudAt9rC0RizKwM9ZcqUfgcddNDw5557rr8Z3TXWWGPuscce+9r++++/5AQKM562GnzVVVet9N57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" width="720">


grouping data
-------------

As you can see, due to the overlapping of different instruments and to
different time snapshots, some points have multiple values. Although
this is not a problem for the fit process, you might want to rebin your
data. This can be obtained with the following command:

.. code:: ipython3

    sed_data.group_data(bin_width=0.2)


.. parsed-literal::

    ===================================================================================================================
    
    ***  binning data  ***
    ---> N bins= 89
    ---> bin_widht= 0.2
    ===================================================================================================================
    


handling errors and systematics
-------------------------------

Another important issues when dealing with fitting of data, is the
proper handling of errors. Typically one might need to add systematics
for different reasons:

-  data are not really simultaneous, and you want to add systematics to
   take this into account
-  data (typically IR up to UV), might have very small errors compared
   to those at higher energies. This might bias the minimizer to
   accomodate the parameters in order to fit 'better' the low
   frequencies branch.

For these reasons the package offer the possibility to add systematics

.. code:: ipython3

    %matplotlib inline
    sed_data.add_systematics(0.2,[10.**6,10.**29])
    myPlot=sed_data.plot_sed()



.. image:: Jet_example_load_data_files/Jet_example_load_data_55_0.png


with this command we add 20% systematics for data between :math:`10^{6}<\nu<10^{29}` Hz
