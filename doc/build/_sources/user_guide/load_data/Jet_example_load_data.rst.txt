
.. _data-format:

Data
====

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
data format is to build an empyt table to have a look at the structure
of the table:

.. code:: ipython3

    from jetset.data_loader import Data
    data=Data(n_rows=10)

we can easily access the astropy table

.. code:: ipython3

    data.table




.. raw:: html

    <i>Table length=10</i>
    <table id="table90625659736" class="table-striped table-bordered table-condensed">
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

    ['/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-stable-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_3C345.dat',
     '/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-stable-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk421.dat',
     '/Users/orion/anaconda3/lib/python3.7/site-packages/jetset-stable-py3.7-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk501.dat']



As you can see there are three 3 files. We use in our example the file for Mrk 421, and we use class:`jetset.data_loader.Data` class  

.. code:: ipython3

    from jetset.data_loader import Data

.. code:: ipython3

    data=Data(data_table=test_SEDs[1])

.. code:: ipython3

    data.table




.. raw:: html

    <i>Table length=112</i>
    <table id="table90625659904" class="table-striped table-bordered table-condensed">
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

    <i>Table length=112</i>
    <table id="table90653823048" class="table-striped table-bordered table-condensed">
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

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAtAAAAGwCAYAAACAS1JbAAAgAElEQVR4Xu2dB7hdRbn+35MeCBBKMBAJQSUJNQmJdIGASg2gNEW6V4pcLEgv3ni9GBC4gghcojTlj5ALKEWQXkTKlUBIaAFpIiDSAoSE9P/zxXVgZ2fvvcqstWbWPr/9PHkIZ89837d+7+zMe2bPmtUhXhCAAAQgAAEIQAACEIBAYgIdiVvSEAIQgAAEIAABCEAAAhAQBppBAAEIQAACEIAABCAAgRQEMNApYNEUAhCAAAQgAAEIQAACGGjGAAQgAAEIQAACEIAABFIQwECngEVTCEAAAhCAAAQgAAEIYKAZAxCAAAQgAAEIQAACEEhBAAOdAhZNIQABCEAAAhCAAAQggIFmDEAAAhCAAAQgAAEIQCAFAQx0Clg0hQAEIAABCEAAAhCAAAaaMQABCEAAAhCAAAQgAIEUBDDQKWDRFAIQgAAEIAABCEAAAhhoxgAEIAABCEAAAhCAAARSEMBAp4BFUwhAAAIQgAAEIAABCGCgGQMQgAAEIAABCEAAAhBIQQADnQIWTSEAAQhAAAIQgAAEIICBZgxAAAIQgAAEIAABCEAgBQEMdApYNIUABCAAAQhAAAIQgAAGmjEAAQhAAAIQgAAEIACBFAQw0Clg0RQCEIAABCAAAQhAAAIYaMYABCAAAQhAAAIQgAAEUhDAQKeARVMIQAACEIAABCAAAQhgoBkDEIBA6QSefPLJXnPnzv1lR0fHlpK6l14ACSGwJIEFixYtur9Xr17fWm+99eYCBwIQgEAcAQx0HCHehwAEcifw2GOPfbd///5HDR48+L1u3botyj0BASGQgsDChQs7Xn755f7vvffez0eNGnVuiq40hQAEuigBDHQXFZ7LhoBPAlOmTJkyfPjwjt69e8/zWQe5IdBJYM6cOT2nT5++cMSIEaOgAgEIQCCOAAY6jhDvQwACuROYMmXKixtuuOHbHR38E5Q7XAJmIrBo0SJNnTp1pZEjR34mUwA6QQACXYoAs1eXkpuLhUAYBKZMmfLSiBEj3gqjGqqAwL8IPP7446uMHDlyCDwgAAEIxBHAQMcR4n0IQCB3Aq4Gerdf3D/Mirr+37ecnntxBCydwM9//vOVL7/88lUmT57sVU8MdOnSkxAClSWAga6sdBQOgeoSCN1Ab7LJJsP22Weft4888si3d9ttt7WmTZu27GuvvdbrxhtvfHaXXXb5oJP8DTfcsNxpp5222lNPPbXM8ssvv+DVV1+dVqvK9OnTex1wwAFDpk6duuzAgQPn/uxnP/vb7rvv/nH/zrabbrrp0Icffni5uXPnTu7Zs2cqYTfZZJOhzz77bN958+Z1GzRo0JxTTz31tf32229GZ5DTTjtt1QsuuOBT7733Xo8hQ4Z8dPbZZ7+y/fbbz7T3W9X/6quv9jj88MPXsLpmz57dbe211/7orLPOemXbbbf9sFGBe+yxx5Drrrtu5SuuuOKv3/jGN97rbHPIIYescemll6567rnnvvSd73zn7UZ9kxroPffcc8i111678rRp055Yf/3158yePbvjwAMPHHz//fcvb9c3ePDgOePHj//7Pvvs877lmTNnTkcr/eprwUCnGno0hkCXJoCB7tLyc/EQ8EOgSgb6jDPOGLDJJpvM2m+//T5z6aWXvlhroO+6665lnn766T5mMH/2s5+tVm+gR44cOXzMmDEzzznnnFevueaaFY466qgh06dPf2L11Vef30n+wgsvXOniiy8eMHny5H7NDLSZ06233vqDRgb04Ycf7rvRRhvNNuN91113LTtu3LihTz311BNrrrnmPPv/nXfeeehtt902fYsttph15plnDjj99NNXf/PNNx/v0aOHtW9a/1NPPdVr0qRJKx588MHvDBo0aN4555yzymmnnTbopZdemrbCCissrB85VqNdw7Bhw2bfeuutz9v78+bN06BBgzbs06fPwmOOOeb1RvVbmwsvvDB2BfrWW2/td/LJJw+yHJ0G+v333+82fvz4gYcddthba6+99txJkyat8M1vfvMzjz766JPDhg2bawa6lX4YaD+ff7JCoB0IYKDbQUWuAQIVI+BioK946OWVfnTjk0PmLVjUsepyved+Z7u1X91v0zXfyRNB5wr00Ucf/fE+7U996lMbXnzxxUsY6M6cv//975c78sgjh9Qa6KlTp/YeM2bMem+88caUFVdccbHhHD16tK1sv3Pccce9af//zjvvdN9oo43WueSSS17cbrvthmcx0LXXbYZ4xx13HH7LLbc8s+22286aOHHiiuedd97AadOmPW3tzHCusMIKo1566aWpZrBb1d+IZ79+/Ubdcsst07/whS/MamSgV1555fm2Cv30008/MWDAgAVXXXXVChdccMGqH374YbcDDzzwLTPQttp82WWXDRg1atSH11xzzcoHHnjgPz/3uc/Nqd3Ccdhhh336kUceWfb222//60orrbTATPaGG2647mWXXfbipptuum6ngW5U49ChQ9c96aSTXjvooIM+XoW3dq3064zDCnSenyJiQaC9CWCg21tfrg4CQRLIaqDNPP/4pqfWnDN/YbfOC+vdo9vCU3dZ9+U8TXQeBvrXv/51//Hjxw964YUXnuys9YADDhjc0dGx6PLLL3/Ffrb//vsP/uxnP/vR17/+9RnDhw/fIKuBHjt27OceeOCB5efOndux5ZZbvn/vvfc+161bN7399tvdttpqq2HnnXfe37baaqsPzzjjjFWvuOKKVZ588smn7P00BvrBBx/su80226zz+uuvP26mtpGBHjRo0Nw333yz58iRI2cdf/zxb+60006f2W233WZMnDhxQK2BPvroo4f86Ec/euW44477p60SX3LJJSuZgf7LX/4yfd99913z73//e6+bb775+eWXX37xLx6nnnrqp+xnl1566SsdHR2jmxnoV155pcfnPve5DR966KGnRo0a9VFtjRjoIP8poCgIVJYABrqy0lE4BKpLII2B7rxh0K72qdffX9ZWnuuvvGf3jkXrrrb84r25edxYmIeBPv/881eaOHHiqo8//vgznfUeddRRg1577bWe11577Uv33XffMoceeuiQadOmPfXCCy/0cjHQFt+MqK2E25aS8ePH/9N+tnDhQp100kkDzzrrrNUldSy33HLzf//73z+39dZbL7GC3GgFvZaxGfHNNtts+B577PHOhAkT/tFo5NkWDjPQO++88/vHH3/8p++6667n1l577fVffvnlqVtvvfXQWgM9YcKE1V9//fWP94vbqvSvfvWrAWusscbc+fPnd9xwww0v9O7de/EDdp577rme22233bApU6Y8bca9mYG26x87duzaQ4YMmXPllVe+XF8jBrq6/15QOQRCJICBDlEVaoJAmxPIaqAf//t7/ZqhGfHpFRbfGBeKgbYV6B/96EeDnn/++Y9XoA888MA1rEZbSR0xYsQ6EyZMeGWXXXaZaTcb1hto24rw+uuv97L2H330UbcePXossj/2/7vvvvs7v/nNb/7WiMUXvvCFtQ8//PB/2o18Z5999irnnnvuwJtuuuk5u+nuuuuuW/6www5ba/LkyU8NGTIk0RaOmTNnmjEdaivlV1111VLGtLOGTgP985///LXBgwevv+OOO8549913u0+aNOll27pSv4Xj0Ucf/fgXCzPQP/zhD9ewveT33Xff05tsssnszrjbb7/9Z3fdddcZdkOn/ayRgbZfFHbdddfPfPDBB91uu+225zvNdy0fDHSb/6PC5UGgZAIY6JKBkw4CEJDSGOhaXhufdscG//xgzmJTWfuyvdD/d/IXlzgBw4VzHivQtgf685///Hr/+Mc/Pt4DPWbMmGF77733O4cccsg7q6666siVVlpp8c2ECxYs0IwZM3rYHuIrrrji+R122GHxLwO15rTZTYT117n55psP/fKXvzzDVqFty0jPnj0XXXzxxYu3jNhr+PDh6x5//PGvH3zwwe92/qzZCrSdcvGlL33pcyuuuOL866+//sXabR/1eWsN9NFHH736Oeecs9oNN9ww3X5BqDfQ9UfWdZ7Cceihh7556qmnfvr222+fPmLEiDmWY7nllhvZq1evRZ0P3Xn77bd79O/ff/5PfvKTV4444oh3zDzvvffeQ1555ZVed95553P9+vVr+Gh4DLTLJ4K+EIBAPQEMNGMCAhAonUBWA+1jD7SZSHtK3VprrbXBhRde+NIOO+zwQZ8+fRaZmTTz9tFHH3XceOONy33/+99f869//esT3bt3V+cK6IgRI4ZvsskmH5/CYTca2ikcAwcOnG/HxHWCty0ctr/4hRdemGondNSvoDY7heOxxx7r8+yzz/baeeedP7BTOC6++OIVv/vd7w65++67n9lyyy1nnXfeeSufffbZq918883PDh8+fO7111+//L777vvZBx544GnbI9yqftsSseOOO362W7dui2655Zbn64/X61w1f+aZZ6bZiRe1BvqNN97o/tBDDy0zbty4D4xTUgNt50BbzXbaxx133DHdVs2Nk9XZ+Ro8ePCIO++885mNN954lpnlfffdd/CTTz65zH333fdso9NBWulXP/C5ibD0fwpICIHKEsBAV1Y6CodAdQlkNdB2xWWdwvG1r33t7e9///tvDRo0aAM7A7qWdqdpvOmmm5azY+Nq39t4441nPvzww4sfCBKdA73W448/vvgc6HPOOeflRudAN9rCURuzmYGePHlyn4MPPnjI888/39eM7pprrjnnuOOOe/2AAw5YfAKFGU9bDb766qtXfv/993t8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" width="720">


**to have interactive plot in ipython**

.. code:: ipython3

    from matplotlib import pylab as plt
    plt.ion()
    myPlot=sed_data.plot_sed()



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAtAAAAGwCAYAAACAS1JbAAAgAElEQVR4Xu2dB7hdRbn+35MeCBBKMBAJQSUJNQmJdIGASg2gNEW6V4pcLEgv3ni9GBC4gghcojTlj5ALKEWQXkTKlUBIaAFpIiDSAoSE9P/zxXVgZ2fvvcqstWbWPr/9PHkIZ89837d+7+zMe2bPmtUhXhCAAAQgAAEIQAACEIBAYgIdiVvSEAIQgAAEIAABCEAAAhAQBppBAAEIQAACEIAABCAAgRQEMNApYNEUAhCAAAQgAAEIQAACGGjGAAQgAAEIQAACEIAABFIQwECngEVTCEAAAhCAAAQgAAEIYKAZAxCAAAQgAAEIQAACEEhBAAOdAhZNIQABCEAAAhCAAAQggIFmDEAAAhCAAAQgAAEIQCAFAQx0Clg0hQAEIAABCEAAAhCAAAaaMQABCEAAAhCAAAQgAIEUBDDQKWDRFAIQgAAEIAABCEAAAhhoxgAEIAABCEAAAhCAAARSEMBAp4BFUwhAAAIQgAAEIAABCGCgGQMQgAAEIAABCEAAAhBIQQADnQIWTSEAAQhAAAIQgAAEIICBZgxAAAIQgAAEIAABCEAgBQEMdApYNIUABCAAAQhAAAIQgAAGmjEAAQhAAAIQgAAEIACBFAQw0Clg0RQCEIAABCAAAQhAAAIYaMYABCAAAQhAAAIQgAAEUhDAQKeARVMIQAACEIAABCAAAQhgoBkDEIBA6QSefPLJXnPnzv1lR0fHlpK6l14ACSGwJIEFixYtur9Xr17fWm+99eYCBwIQgEAcAQx0HCHehwAEcifw2GOPfbd///5HDR48+L1u3botyj0BASGQgsDChQs7Xn755f7vvffez0eNGnVuiq40hQAEuigBDHQXFZ7LhoBPAlOmTJkyfPjwjt69e8/zWQe5IdBJYM6cOT2nT5++cMSIEaOgAgEIQCCOAAY6jhDvQwACuROYMmXKixtuuOHbHR38E5Q7XAJmIrBo0SJNnTp1pZEjR34mUwA6QQACXYoAs1eXkpuLhUAYBKZMmfLSiBEj3gqjGqqAwL8IPP7446uMHDlyCDwgAAEIxBHAQMcR4n0IQCB3Aq4Gerdf3D/Mirr+37ecnntxBCydwM9//vOVL7/88lUmT57sVU8MdOnSkxAClSWAga6sdBQOgeoSCN1Ab7LJJsP22Weft4888si3d9ttt7WmTZu27GuvvdbrxhtvfHaXXXb5oJP8DTfcsNxpp5222lNPPbXM8ssvv+DVV1+dVqvK9OnTex1wwAFDpk6duuzAgQPn/uxnP/vb7rvv/nH/zrabbrrp0Icffni5uXPnTu7Zs2cqYTfZZJOhzz77bN958+Z1GzRo0JxTTz31tf32229GZ5DTTjtt1QsuuOBT7733Xo8hQ4Z8dPbZZ7+y/fbbz7T3W9X/6quv9jj88MPXsLpmz57dbe211/7orLPOemXbbbf9sFGBe+yxx5Drrrtu5SuuuOKv3/jGN97rbHPIIYescemll6567rnnvvSd73zn7UZ9kxroPffcc8i111678rRp055Yf/3158yePbvjwAMPHHz//fcvb9c3ePDgOePHj//7Pvvs877lmTNnTkcr/eprwUCnGno0hkCXJoCB7tLyc/EQ8EOgSgb6jDPOGLDJJpvM2m+//T5z6aWXvlhroO+6665lnn766T5mMH/2s5+tVm+gR44cOXzMmDEzzznnnFevueaaFY466qgh06dPf2L11Vef30n+wgsvXOniiy8eMHny5H7NDLSZ06233vqDRgb04Ycf7rvRRhvNNuN91113LTtu3LihTz311BNrrrnmPPv/nXfeeehtt902fYsttph15plnDjj99NNXf/PNNx/v0aOHtW9a/1NPPdVr0qRJKx588MHvDBo0aN4555yzymmnnTbopZdemrbCCissrB85VqNdw7Bhw2bfeuutz9v78+bN06BBgzbs06fPwmOOOeb1RvVbmwsvvDB2BfrWW2/td/LJJw+yHJ0G+v333+82fvz4gYcddthba6+99txJkyat8M1vfvMzjz766JPDhg2bawa6lX4YaD+ff7JCoB0IYKDbQUWuAQIVI+BioK946OWVfnTjk0PmLVjUsepyved+Z7u1X91v0zXfyRNB5wr00Ucf/fE+7U996lMbXnzxxUsY6M6cv//975c78sgjh9Qa6KlTp/YeM2bMem+88caUFVdccbHhHD16tK1sv3Pccce9af//zjvvdN9oo43WueSSS17cbrvthmcx0LXXbYZ4xx13HH7LLbc8s+22286aOHHiiuedd97AadOmPW3tzHCusMIKo1566aWpZrBb1d+IZ79+/Ubdcsst07/whS/MamSgV1555fm2Cv30008/MWDAgAVXXXXVChdccMGqH374YbcDDzzwLTPQttp82WWXDRg1atSH11xzzcoHHnjgPz/3uc/Nqd3Ccdhhh336kUceWfb222//60orrbTATPaGG2647mWXXfbipptuum6ngW5U49ChQ9c96aSTXjvooIM+XoW3dq3064zDCnSenyJiQaC9CWCg21tfrg4CQRLIaqDNPP/4pqfWnDN/YbfOC+vdo9vCU3dZ9+U8TXQeBvrXv/51//Hjxw964YUXnuys9YADDhjc0dGx6PLLL3/Ffrb//vsP/uxnP/vR17/+9RnDhw/fIKuBHjt27OceeOCB5efOndux5ZZbvn/vvfc+161bN7399tvdttpqq2HnnXfe37baaqsPzzjjjFWvuOKKVZ588smn7P00BvrBBx/su80226zz+uuvP26mtpGBHjRo0Nw333yz58iRI2cdf/zxb+60006f2W233WZMnDhxQK2BPvroo4f86Ec/euW44477p60SX3LJJSuZgf7LX/4yfd99913z73//e6+bb775+eWXX37xLx6nnnrqp+xnl1566SsdHR2jmxnoV155pcfnPve5DR966KGnRo0a9VFtjRjoIP8poCgIVJYABrqy0lE4BKpLII2B7rxh0K72qdffX9ZWnuuvvGf3jkXrrrb84r25edxYmIeBPv/881eaOHHiqo8//vgznfUeddRRg1577bWe11577Uv33XffMoceeuiQadOmPfXCCy/0cjHQFt+MqK2E25aS8ePH/9N+tnDhQp100kkDzzrrrNUldSy33HLzf//73z+39dZbL7GC3GgFvZaxGfHNNtts+B577PHOhAkT/tFo5NkWDjPQO++88/vHH3/8p++6667n1l577fVffvnlqVtvvfXQWgM9YcKE1V9//fWP94vbqvSvfvWrAWusscbc+fPnd9xwww0v9O7de/EDdp577rme22233bApU6Y8bca9mYG26x87duzaQ4YMmXPllVe+XF8jBrq6/15QOQRCJICBDlEVaoJAmxPIaqAf//t7/ZqhGfHpFRbfGBeKgbYV6B/96EeDnn/++Y9XoA888MA1rEZbSR0xYsQ6EyZMeGWXXXaZaTcb1hto24rw+uuv97L2H330UbcePXossj/2/7vvvvs7v/nNb/7WiMUXvvCFtQ8//PB/2o18Z5999irnnnvuwJtuuuk5u+nuuuuuW/6www5ba/LkyU8NGTIk0RaOmTNnmjEdaivlV1111VLGtLOGTgP985///LXBgwevv+OOO8549913u0+aNOll27pSv4Xj0Ucf/fgXCzPQP/zhD9ewveT33Xff05tsssnszrjbb7/9Z3fdddcZdkOn/ayRgbZfFHbdddfPfPDBB91uu+225zvNdy0fDHSb/6PC5UGgZAIY6JKBkw4CEJDSGOhaXhufdscG//xgzmJTWfuyvdD/d/IXlzgBw4VzHivQtgf685///Hr/+Mc/Pt4DPWbMmGF77733O4cccsg7q6666siVVlpp8c2ECxYs0IwZM3rYHuIrrrji+R122GHxLwO15rTZTYT117n55psP/fKXvzzDVqFty0jPnj0XXXzxxYu3jNhr+PDh6x5//PGvH3zwwe92/qzZCrSdcvGlL33pcyuuuOL866+//sXabR/1eWsN9NFHH736Oeecs9oNN9ww3X5BqDfQ9UfWdZ7Cceihh7556qmnfvr222+fPmLEiDmWY7nllhvZq1evRZ0P3Xn77bd79O/ff/5PfvKTV4444oh3zDzvvffeQ1555ZVed95553P9+vVr+Gh4DLTLJ4K+EIBAPQEMNGMCAhAonUBWA+1jD7SZSHtK3VprrbXBhRde+NIOO+zwQZ8+fRaZmTTz9tFHH3XceOONy33/+99f869//esT3bt3V+cK6IgRI4ZvsskmH5/CYTca2ikcAwcOnG/HxHWCty0ctr/4hRdemGondNSvoDY7heOxxx7r8+yzz/baeeedP7BTOC6++OIVv/vd7w65++67n9lyyy1nnXfeeSufffbZq918883PDh8+fO7111+//L777vvZBx544GnbI9yqftsSseOOO362W7dui2655Zbn64/X61w1f+aZZ6bZiRe1BvqNN97o/tBDDy0zbty4D4xTUgNt50BbzXbaxx133DHdVs2Nk9XZ+Ro8ePCIO++885mNN954lpnlfffdd/CTTz65zH333fdso9NBWulXP/C5ibD0fwpICIHKEsBAV1Y6CodAdQlkNdB2xWWdwvG1r33t7e9///tvDRo0aAM7A7qWdqdpvOmmm5azY+Nq39t4441nPvzww4sfCBKdA73W448/vvgc6HPOOeflRudAN9rCURuzmYGePHlyn4MPPnjI888/39eM7pprrjnnuOOOe/2AAw5YfAKFGU9bDb766qtXfv/993t8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" width="720">


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
