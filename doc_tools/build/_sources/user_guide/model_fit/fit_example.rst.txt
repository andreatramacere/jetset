
.. toctree::


.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')
    
    import matplotlib.pylab as plt
    %matplotlib inline

Fit User guide
==============

data format
-----------

The SED data can be stored in ASCII file in a quite flexible way, but
some requirements are needed:

-  you must provide at least two columns for frequencies and fluxes
-  frequencies are in *Hz*
-  fluxes are in *cgs*, .....

The header of the file can contain some meta-data that are sourced when
the data are loaded. The meta-data available are :

-  z : redhsift
-  resframe: restframe of the data ``src`` or ``obs``
-  data\_scale: scale of the data ``lin-lin`` or ``log-log``
-  dataType: structure of the comumns with the SED data

the meaning of these meta-data is explained in detail in :class:`jetset.data_loader.ObsData` class 
documentation. The meta-data can be included in the header with  line like:


.. code:: ipython3

    # metadata
    # md z  0.0308
    # md restframe  obs
    # md data_scale  lin-lin
    # md col_types x,y,dy,data_set
    # md obj_name     J1104+3812,Mrk421
    #

A typical structure of SED data file, including meta-data declaration is
the following:

::

    # metadata
    # md z  0.0308
    # md restframe  obs
    # md data_scale  lin-lin
    # md col_types x,y,dy,data_set
    # md obj_name     J1104+3812,Mrk421
    #
    # Frequency [Hz]  EnergyFlux [erg/cm2/s]
    #  Xval           Yval       YvalError    data_set-flag
    2.299540e+09 1.340900e-14 3.910000e-16    campaing-2009
    2.639697e+09 1.793088e-14 3.231099e-26    campaing-2009
    4.799040e+09 2.313600e-14 2.400000e-16    campaing-2009
    4.805039e+09 1.773414e-14 1.773414e-15    campaing-2009
    4.843552e+09 2.776140e-14 2.615339e-26    campaing-2009
    7.698460e+09 3.696000e-14 4.620000e-16    campaing-2009
    8.267346e+09 2.836267e-14 2.836267e-15    campaing-2009
    8.331867e+09 3.989630e-14 3.627671e-26    campaing-2009
    8.388659e+09 3.163450e-14 1.931495e-15    campaing-2009
    8.399994e+09 4.000500e-14 5.041094e-15    campaing-2009
    1.044892e+10 4.626737e-14 3.297726e-26    campaing-2009
    1.109778e+10 4.617600e-14 6.660000e-16    campaing-2009
    1.456571e+10 5.628417e-14 4.453463e-26    campaing-2009

Loading SEDs
------------


The most effective way to import the SED data is to create an object 
instance of :class:`jetset.data_loader.ObsData` class 
(see the documentation for the :doc:`data_loader <../../modules_doc/data_loader>` module)
The package provides some test SEDs, accessible as follows:

.. code:: ipython3

    from jetset.test_data_helper import  test_SEDs
    test_SEDs




.. parsed-literal::

    ['/Users/orion/anaconda2/envs/py36/lib/python3.6/site-packages/jetset-1.2.0-py3.6-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_3C345.dat',
     '/Users/orion/anaconda2/envs/py36/lib/python3.6/site-packages/jetset-1.2.0-py3.6-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk421.dat',
     '/Users/orion/anaconda2/envs/py36/lib/python3.6/site-packages/jetset-1.2.0-py3.6-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk501.dat']



to load the SED of Mrk 421, the first one in the list:

.. code:: ipython3

    mySED=test_SEDs[1]
    from jetset.data_loader import ObsData
    sed_data=ObsData(data_file=mySED)


.. parsed-literal::

    =============================================================================================
    
    *** getting meta-data from file header
    col_types None
    set md z  to 0.0308
    set md restframe  to obs
    set md data_scale  to lin-lin
    set md col_types  to x,y,dy,data_set
    set md obj_name  to J1104+3812,Mrk421
    ciccio x,y,dy,data_set
    col_types x,y,dy,data_set
    =============================================================================================
    
    col_types a x,y,dy,data_set
    col_types c x,y,dy,data_set
    col_types b ['x', 'y', 'dy', 'data_set']
    ciccio [('nu_data', 'f8'), ('dnu_data', 'f8'), ('nuFnu_data', 'f8'), ('dnuFnu_data', 'f8'), ('nu_data_log', 'f8'), ('dnu_data_log', 'f8'), ('nuFnu_data_log', 'f8'), ('dnuFnu_data_log', 'f8'), ('dnuFnu_facke', 'f8'), ('dnuFnu_facke_log', 'f8'), ('UL', 'bool'), ('zero_error', 'bool'), ('T_start', 'f8'), ('T_stop', 'f8'), ('data_set', 'S16')]
    =============================================================================================
    
    *** loading data ***
    ---> loading data for file=/Users/orion/anaconda2/envs/py36/lib/python3.6/site-packages/jetset-1.2.0-py3.6-macosx-10.7-x86_64.egg/jetset/test_data/SEDs_data/SED_MW_Mrk421.dat
    ---> found these col ID=range(0, 4) and names=['x', 'y', 'dy', 'data_set']:
    ---> z=3.080000e-02
    ---> restframe=obs
    ---> obj_name=J1104+3812,Mrk421 
    ---> data_scale=lin-lin 
    col_types ['x', 'y', 'dy', 'data_set']
    ['x', 'y', 'dy', 'data_set'] range(0, 4)
    nu_data x
    nuFnu_data y
    dnuFnu_data dy
    data_set data_set
    ---> data len=112
    ---> Settin  UL for val 0
    ---> Settin  UL for val 0.2
    =============================================================================================
    


As you can see the all the meta-data have been properly sourced from the
SED file header. You also get information on the lenght of the data,
before and after elimination of duplicated entries, and upper limits
These meta-data are parameters needed by the

:class:`jetset.data_loader.ObsData` constructor. 

Plotting data
-------------

We can now plot our SED using the :class:`BlazarSEDFit.plot_sedfit.Plot` class 
(see the documentation for the :doc:`plot_sedfit <../../modules_doc/plot_sedfit>` module)

.. code:: ipython3

    from jetset.plot_sedfit import Plot
    myPlot=Plot(sed_data,interactive=True)
    
    myPlot.add_data_plot(sed_data,autoscale=True)


.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing



.. image:: fit_example_files/fit_example_17_1.png


grouping data
-------------

As you can see, due to the overlapping of different instruments and to
different time snapshots, some points have multiple values. Although
this is not a problem for the fit process, you might want to rebin your
data. This can be obtained with the following command:

.. code:: ipython3

    sed_data.group_data(bin_width=0.2)


.. parsed-literal::

    =============================================================================================
    
    ***  binning data  ***
    ---> N bins= 89
    ---> bin_widht= 0.2
    =============================================================================================
    


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

    sed_data.add_systematics(0.2,[10.**6,10.**29])
    myPlot=Plot(interactive=True)
    myPlot.add_data_plot(sed_data,label='grouped+syst',autoscale=True)


.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing



.. image:: fit_example_files/fit_example_21_1.png


with this command we add 20% systematics for data between :math:`10^{6}<\nu<10^{29}` Hz

SEDShape: Spectral indices
--------------------------

.. code:: ipython3

    from jetset.sed_shaper import  SEDShape
    my_shape=SEDShape(sed_data)
    my_shape.eval_indices()


.. parsed-literal::

    =============================================================================================
    
    *** evaluating spectral indices for data ***
    ---> range for indexradio updated to [6.000000,10.000000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=80, res=0.040981, chisq=0.36841**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 0.36841838886162864 | total call = 87 | ncalls = 87
    edm = 1.8023341884611675e-07 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.5973  |  0.2061  |          |          | -10      |  10      |          |
    |    1 | par_1 = -19.36   |  1.992   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.5973317132652589, -19.35835367872672)
    res check 0.0380162913807775 0.3684200285851512
    ---> 1000000.0 10000000000.0 100
    ---> name = radio            range=[6.000 ,10.000] log(Hz)  photon.val=-1.402671e+00, err=2.061105e-01 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +5.973286e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.935838e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=88
    mesg=
    dof=1
    chisq=0.368420, chisq/red=0.368420 null hypothesis sig=0.543867
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +5.973286e-01 | +2.061105e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.935838e+01 | +1.991979e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexradio_mm updated to [10.000000,11.000000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=120, res=-0.000271, chisq=0.01395**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 0.013955396221108125 | total call = 124 | ncalls = 124
    edm = 9.167821380589557e-06 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.7106  |  0.3008  |          |          | -10      |  10      |          |
    |    1 | par_1 = -20.49   |  3.102   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.710563240761811, -20.48609688874422)
    res check -0.0018041570213917016 0.013956525602712383
    ---> 10000000000.0 100000000000.0 100
    ---> name = radio_mm         range=[10.000,11.000] log(Hz)  photon.val=-1.289439e+00, err=3.008381e-01 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +7.105609e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -2.048612e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=125
    mesg=
    dof=1
    chisq=0.013957, chisq/red=0.013957 null hypothesis sig=0.905959
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +7.105609e-01 | +3.008381e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -2.048612e+01 | +3.101753e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexmm_IR updated to [10.300000,13.700000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=100, res=0.002166, chisq=0.16989**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 0.1698919441697112 | total call = 108 | ncalls = 108
    edm = 2.4889897383547027e-06 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.8903  |  0.1269  |          |          | -10      |  10      |          |
    |    1 | par_1 = -22.37   |  1.363   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.8902738412904547, -22.367394988333345)
    res check 0.0012528014355625283 0.16989620992207333
    ---> 19952623149.68891 50118723362726.945 100
    ---> name = mm_IR            range=[10.300,13.700] log(Hz)  photon.val=-1.109729e+00, err=1.268819e-01 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +8.902713e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -2.236742e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=109
    mesg=
    dof=1
    chisq=0.169896, chisq/red=0.169896 null hypothesis sig=0.680204
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +8.902713e-01 | +1.268819e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -2.236742e+01 | +1.363095e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexIR_Opt updated to [12.500000,14.500000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=150, res=-0.063099, chisq=0.03255**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 0.03255413239842313 | total call = 155 | ncalls = 155
    edm = 5.5388584543726864e-05 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.2179  |  0.4602  |          |          | -10      |  10      |          |
    |    1 | par_1 = -13.33   |  6.362   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.21786160571751623, -13.328291691141612)
    res check -0.0656184764763343 0.032552626185190045
    ---> 3162277660168.392 316227766016836.6 100
    ---> name = IR_Opt           range=[12.500,14.500] log(Hz)  photon.val=-1.782141e+00, err=4.602299e-01 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +2.178594e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.332832e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=156
    mesg=
    dof=1
    chisq=0.032553, chisq/red=0.032553 null hypothesis sig=0.856820
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +2.178594e-01 | +4.602299e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.332832e+01 | +6.361702e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexOpt_UV updated to [14.000000,16.000000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=110, res=-0.044237, chisq=1.17052**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 1.1705275574077267 | total call = 119 | ncalls = 119
    edm = 7.233067985143407e-09 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.3788  |  0.09841 |          |          | -10      |  10      |          |
    |    1 | par_1 = -15.64   |  1.449   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.37882181023299566, -15.636958496747097)
    minim function calls=120, res=-0.047985, chisq=1.17052res check -0.047984938876328154 1.170529382549224
    ---> 100000000000000.0 1e+16 100
    ---> name = Opt_UV           range=[14.000,16.000] log(Hz)  photon.val=-1.621180e+00, err=9.841047e-02 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +3.788199e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.563699e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=120
    mesg=
    dof=5
    chisq=1.170529, chisq/red=0.234106 null hypothesis sig=0.947673
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +3.788199e-01 | +9.841047e-02 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.563699e+01 | +1.449001e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexBBB updated to [14.800000,16.200000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=200, res=0.023204, chisq=0.16084**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 0.16083664912879314 | total call = 200 | ncalls = 200
    edm = 1.3739884290102432e-05 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.7254  |  0.3847  |          |          | -10      |  10      |          |
    |    1 | par_1 = -20.86   |  5.63    |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.7253961177267225, -20.86309161951963)
    res check 0.023203710060429455 0.16084095856515768
    ---> 630957344480194.2 1.5848931924611238e+16 100
    ---> name = BBB              range=[14.800,16.200] log(Hz)  photon.val=-1.274606e+00, err=3.846946e-01 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +7.253943e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -2.086312e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=201
    mesg=
    dof=1
    chisq=0.160841, chisq/red=0.160841 null hypothesis sig=0.688383
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +7.253943e-01 | +3.846946e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -2.086312e+01 | +5.629936e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexUV_X updated to [15.000000,17.500000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=90, res=-0.013896, chisq=1.13058**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 1.1305823072504297 | total call = 91 | ncalls = 91
    edm = 7.65798990257511e-06 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.1541  |  0.03713 |          |          | -10      |  10      |          |
    |    1 | par_1 = -12.19   |  0.6164  |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.15408070210057012, -12.187409638022114)
    res check -0.017728961580666613 1.1305788936641659
    ---> 1000000000000000.0 3.1622776601683795e+17 100
    ---> name = UV_X             range=[15.000,17.500] log(Hz)  photon.val=-1.845921e+00, err=3.712945e-02 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +1.540792e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.218744e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=92
    mesg=
    dof=4
    chisq=1.130579, chisq/red=0.282645 null hypothesis sig=0.889391
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +1.540792e-01 | +3.712945e-02 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.218744e+01 | +6.163737e-01 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexX updated to [16.000000,19.000000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=110, res=-1.103361, chisq=21.93071**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 21.93070382034992 | total call = 117 | ncalls = 117
    edm = 6.068854602007761e-05 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 = -0.4578  |  0.04767 |          |          | -10      |  10      |          |
    |    1 | par_1 = -1.606   |  0.8467  |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (-0.45779697386609364, -1.606322923127852)
    res check -1.126465635005074 21.930709522488655
    ---> 1e+16 1e+19 100
    ---> name = X                range=[16.000,19.000] log(Hz)  photon.val=-2.457801e+00, err=4.767152e-02 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | -4.578008e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.606391e+00 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=118
    mesg=
    dof=9
    chisq=21.930710, chisq/red=2.436746 null hypothesis sig=0.009101
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | -4.578008e-01 | +4.767152e-02 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.606391e+00 | +8.467033e-01 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    ---> range for indexFermi updated to [22.380000,25.380000]
    directory ./ already existing
    directory .//spectral-indices-best-fit/ already existing
    removing existing dir
    the directory .//spectral-indices-best-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10.0, 10.0), 'par_1': -10.0, 'limit_par_1': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=110, res=-0.149538, chisq=1.21422**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 1.2142224037471625 | total call = 115 | ncalls = 115
    edm = 2.013143960876923e-06 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.2044  |  0.04438 |          |          | -10      |  10      |          |
    |    1 | par_1 = -15.27   |  1.056   |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.20438294397813905, -15.272890189690562)
    res check -0.1556421183794186 1.214224314105566
    ---> 2.3988329190194848e+22 2.398832919019485e+25 100
    ---> name = Fermi            range=[22.380,25.380] log(Hz)  photon.val=-1.795618e+00, err=4.437842e-02 
    
    **************************************************************************************************
    Fit report
    
    Model: spectral-indices-best-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     alpha            | spectral-slope           |                  | +2.043818e-01 | [-1.000000e+01,+1.000000e+01]  
     K                | flux-const               | erg cm^-2 s^-1   | -1.527292e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=116
    mesg=
    dof=6
    chisq=1.214224, chisq/red=0.202371 null hypothesis sig=0.976176
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     alpha            | +2.043818e-01 | +4.437842e-02 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     K                | -1.527292e+01 | +1.055578e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    
    
    =============================================================================================
    


.. code:: ipython3

    myPlot=Plot(sed_data,interactive=True)
    
    for model in my_shape.index_models:
        myPlot.add_model_plot(model,label=model.name,line_style='--',autoscale=True)
    
    myPlot.add_data_plot(sed_data,autoscale=True,label='data',color='red')
    myPlot.rescale(y_min=-14,y_max=-8,x_min=8)


.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing
    label radio
    label radio_mm
    label mm_IR
    label IR_Opt
    label Opt_UV
    label BBB
    label UV_X
    label X
    label Fermi



.. image:: fit_example_files/fit_example_25_1.png


SEDShape: Log-Log Polynomila fit
--------------------------------

.. code:: ipython3

    myPlot.save('SED_indices_rebinned.png')
    
    
    my_shape.sync_fit(check_host_gal_template=True)


.. parsed-literal::

    =============================================================================================
    
    *** Log-Polynomial fitting of the synchrotron component ***
    ---> first blind fit run, log-cubic fit range: [9, 19]
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     b                | curvature                |                  | -1.000000e+00 | [-1.000000e+01,+0.000000e+00]  
     c                | third-degree             |                  | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     Ep               | peak freq                | Hz               | +1.400000e+01 | [+0.000000e+00,+3.000000e+01]  
     Sp               | peak flux                | erg cm^-2 s^-1   | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    directory ./ already existing
    directory .//sync-shape-fit/ already existing
    removing existing dir
    the directory .//sync-shape-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10, 0), 'par_1': -1.0, 'limit_par_1': (-10.0, 10.0), 'par_2': 14.0, 'limit_par_2': (0.0, 30.0), 'par_3': -10.0, 'limit_par_3': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=800, res=-0.657228, chisq=13.73314**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 13.73314818997152 | total call = 808 | ncalls = 808
    edm = 5.038846082587436e-05 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 = -0.1595  |  0.009421 |          |          | -10      |  0       |          |
    |    1 | par_1 = -0.0109  |  0.001313 |          |          | -10      |  10      |          |
    |    2 | par_2 =  16.7    |  0.04841 |          |          |  0       |  30      |          |
    |    3 | par_3 = -9.485   |  0.03388 |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (-0.15952729965761314, -0.010904800855310981, 16.695980883689245, -9.48449545453073)
    res check -0.6634820658680649 13.733155030005754
    ---> 1000000000.0 1e+19 100
    directory ./ already existing
    directory .//sync-shape-fit/ already existing
    removing existing dir
    the directory .//sync-shape-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10, 0), 'par_1': -1.0, 'limit_par_1': (-10.0, 10.0), 'par_2': 14.0, 'limit_par_2': (0.0, 30.0), 'par_3': -10.0, 'limit_par_3': (-30.0, 0.0), 'par_4': -9.484531267217534, 'limit_par_4': (-11.484531267217534, -7.484531267217534), 'par_5': 0, 'limit_par_5': (-0.5, 0.5)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=1480, res=-0.763216, chisq=13.28537**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 13.285367168384496 | total call = 1488 | ncalls = 1480
    edm = 0.009548166333585998 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |           True |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 = -0.1604  |  2.1     |          |          | -10      |  0       |          |
    |    1 | par_1 = -0.01109 |  0.1674  |          |          | -10      |  10      |          |
    |    2 | par_2 =  16.71   |  9.703   |          |          |  0       |  30      |          |
    |    3 | par_3 = -9.491   |  9.277   |          |          | -30      |  0       |          |
    |    4 | par_4 = -11.07   |  2.097   |          |          | -11.48   | -7.485   |          |
    |    5 | par_5 = -0.01811 |  0.08444 |          |          | -0.5     |  0.5     |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (-0.16036897465601996, -0.0110925794540222, 16.714910037550833, -9.490914998591467, -11.065221264625162, -0.018113910210059703)
    res check -0.7632471345640286 13.285367169522232
    ---> 1000000000.0 1e+19 100
    
    **************************************************************************************************
    Fit report
    
    Model: sync-shape-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     b                | curvature                |                  | -1.603690e-01 | [-1.000000e+01,+0.000000e+00]  
     c                | third-degree             |                  | -1.109258e-02 | [-1.000000e+01,+1.000000e+01]  
     Ep               | peak freq                | Hz               | +1.671491e+01 | [+0.000000e+00,+3.000000e+01]  
     Sp               | peak flux                | erg cm^-2 s^-1   | -9.490915e+00 | [-3.000000e+01,+0.000000e+00]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | -1.811391e-02 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=1489
    mesg=
    dof=19
    chisq=13.285367, chisq/red=0.699230 null hypothesis sig=0.823647
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     b                | -1.603690e-01 | +2.100083e+00 | -1.000000e+00 | [-1.000000e+01,+0.000000e+00]
     c                | -1.109258e-02 | +1.674351e-01 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     Ep               | +1.671491e+01 | +9.702709e+00 | +1.400000e+01 | [+0.000000e+00,+3.000000e+01]
     Sp               | -9.490915e+00 | +9.277498e+00 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
     nuFnu_p_host     | -1.106522e+01 | +2.097009e+00 | -9.484531e+00 | [-1.148453e+01,-7.484531e+00]
     nu_scale         | -1.811391e-02 | +8.443645e-02 | +0.000000e+00 | [-5.000000e-01,+5.000000e-01]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    ---> class:  HSP
    ---> sync       nu_p=+1.671491e+01 (err=+9.702709e+00)  nuFnu_p=-9.490915e+00 (err=+9.277498e+00) curv.=-1.603690e-01 (err=+2.100083e+00)


.. code:: ipython3

    my_shape.IC_fit()


.. parsed-literal::

    =============================================================================================
    
    *** Log-Polynomial fitting of the IC component ***
    ---> log-cubic fit range: [22, 28]
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     b                | curvature                |                  | -1.000000e+00 | [-1.000000e+01,+0.000000e+00]  
     c                | third-degree             |                  | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     Ep               | peak freq                | Hz               | +2.525747e+01 | [+0.000000e+00,+3.000000e+01]  
     Sp               | peak flux                | erg cm^-2 s^-1   | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    directory ./ already existing
    directory .//IC-shape-fit/ already existing
    removing existing dir
    the directory .//IC-shape-fit/ has been created
    dict {'par_0': -1.0, 'limit_par_0': (-10, 0), 'par_1': -1.0, 'limit_par_1': (-10.0, 10.0), 'par_2': 25.25746738962698, 'limit_par_2': (0.0, 30.0), 'par_3': -10.0, 'limit_par_3': (-30.0, 0.0)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=210, res=-0.728174, chisq=3.43994**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 3.4399459548662534 | total call = 218 | ncalls = 218
    edm = 3.283406576818471e-05 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |           True |           True |           True |           True |          False |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 = -0.2056  |  0.04002 |          |          | -10      |  0       |          |
    |    1 | par_1 = -0.04991 |  0.01563 |          |          | -10      |  10      |          |
    |    2 | par_2 =  25.26   |  0.1075  |          |          |  0       |  30      |          |
    |    3 | par_3 = -10.12   |  0.04768 |          |          | -30      |  0       |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (-0.205604968378843, -0.04990591734117622, 25.25922273964055, -10.122349495222544)
    res check -0.7341759774539718 3.439939831614986
    ---> 1e+22 1e+28 100
    
    **************************************************************************************************
    Fit report
    
    Model: IC-shape-fit
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     b                | curvature                |                  | -2.056050e-01 | [-1.000000e+01,+0.000000e+00]  
     c                | third-degree             |                  | -4.990592e-02 | [-1.000000e+01,+1.000000e+01]  
     Ep               | peak freq                | Hz               | +2.525916e+01 | [+0.000000e+00,+3.000000e+01]  
     Sp               | peak flux                | erg cm^-2 s^-1   | -1.012238e+01 | [-3.000000e+01,+0.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=219
    mesg=
    dof=12
    chisq=3.439940, chisq/red=0.286662 null hypothesis sig=0.991560
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     b                | -2.056050e-01 | +4.001839e-02 | -1.000000e+00 | [-1.000000e+01,+0.000000e+00]
     c                | -4.990592e-02 | +1.563202e-02 | -1.000000e+00 | [-1.000000e+01,+1.000000e+01]
     Ep               | +2.525916e+01 | +1.075364e-01 | +2.525747e+01 | [+0.000000e+00,+3.000000e+01]
     Sp               | -1.012238e+01 | +4.767663e-02 | -1.000000e+01 | [-3.000000e+01,+0.000000e+00]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    ---> IC         nu_p=+2.525916e+01 (err=+1.075364e-01)  nuFnu_p=-1.012238e+01 (err=+4.767663e-02) curv.=-2.056050e-01 (err=+4.001839e-02)
    =============================================================================================
    


.. code:: ipython3

    myPlot=Plot(sed_data,interactive=True)
    
    
    myPlot.add_model_plot(my_shape.sync_fit,label='sync, poly-fit')
    
    myPlot.add_model_plot(my_shape.host_gal,label='host-gal')
    
    myPlot.add_model_plot(my_shape.sync_fit_model,label='sync+host, poly-fit')
    myPlot.add_model_plot(my_shape.IC_fit_model,label='IC, poly-fit')
    myPlot.add_data_plot(sed_data,autoscale=True)
    myPlot.rescale(y_min=-14,y_max=-8,x_min=8)


.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing
    <bound method SEDShape.sync_fit of <jetset.sed_shaper.SEDShape object at 0x1511e81dd8>> !!! Error has no SED instance or something wrong in get_model_points()
    label host-gal
    label sync+host, poly-fit
    label IC, poly-fit



.. image:: fit_example_files/fit_example_29_1.png


.. code:: ipython3

    my_shape.show_values()


.. parsed-literal::

    =============================================================================================
    
    *** SEDShape values ***
    ---> spectral inidces values
    ---> name = radio            range=[6.000 ,10.000] log(Hz)  photon.val=-1.402671e+00, err=2.061105e-01 
    ---> name = radio_mm         range=[10.000,11.000] log(Hz)  photon.val=-1.289439e+00, err=3.008381e-01 
    ---> name = mm_IR            range=[10.300,13.700] log(Hz)  photon.val=-1.109729e+00, err=1.268819e-01 
    ---> name = IR_Opt           range=[12.500,14.500] log(Hz)  photon.val=-1.782141e+00, err=4.602299e-01 
    ---> name = Opt_UV           range=[14.000,16.000] log(Hz)  photon.val=-1.621180e+00, err=9.841047e-02 
    ---> name = BBB              range=[14.800,16.200] log(Hz)  photon.val=-1.274606e+00, err=3.846946e-01 
    ---> name = UV_X             range=[15.000,17.500] log(Hz)  photon.val=-1.845921e+00, err=3.712945e-02 
    ---> name = X                range=[16.000,19.000] log(Hz)  photon.val=-2.457801e+00, err=4.767152e-02 
    ---> name = Fermi            range=[22.380,25.380] log(Hz)  photon.val=-1.795618e+00, err=4.437842e-02 
    ---> S/IC peak values
    ---> sync       nu_p=+1.671491e+01 (err=+9.702709e+00)  nuFnu_p=-9.490915e+00 (err=+9.277498e+00) curv.=-1.603690e-01 (err=+2.100083e+00)
    
    ---> IC         nu_p=+2.525916e+01 (err=+1.075364e-01)  nuFnu_p=-1.012238e+01 (err=+4.767663e-02) curv.=-2.056050e-01 (err=+4.001839e-02)
    
    
    =============================================================================================
    


Constraining SSC/EC model
-------------------------

.. code:: ipython3

    from jetset.obs_constrain import ObsConstrain
    
    sed_obspar=ObsConstrain(beaming=25,B_range=[0.01,0.1],distr_e='lppl',t_var_sec=3*86400,nu_cut_IR=9.0E12,SEDShape=my_shape)


.. parsed-literal::

    directory ./ already existing


.. code:: ipython3

    jet_model=sed_obspar.constrain_SSC_model()


.. parsed-literal::

    =============================================================================================
    
    ***  constrains parameters from observable ***
    directory .//obs_constrain_lppl/ already existing
    removing existing dir
    the directory .//obs_constrain_lppl/ has been created
    directory ./ already existing
    directory .//lppl_jet_prod/ already existing
    removing existing dir
    the directory .//lppl_jet_prod/ has been created
    directory .//obs_constrain_lppl/ already existing
    removing existing dir
    the directory .//obs_constrain_lppl/ has been created
    -----------------------------------------------------------------------------------------
    model parameters for jet model:
    electron distribution type = lppl  
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +1.000000e+02 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +2.000000e+00 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +1.000000e+08 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.000000e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +4.000000e-01 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.000000e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +3.000000e+15 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +1.000000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +1.000000e-01 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------
    ---> ***  emitting region parameters  ***
    ---> name = beam_obj          type = beaming                   units =                   val = +2.500000e+01  phys-bounds = [+1.000000e+00,No           ]  
    ---> setting par type redshift, corresponding to par z_cosm
    --->  name = z_cosm            type = redshift                  units =                   val = +3.080000e-02  phys-bounds = [+0.000000e+00,No           ]  
    ---> setting par type magnetic_field, corresponding to par B
    --->  name = B                 type = magnetic_field            units = G                 val = +5.500000e-02  phys-bounds = [+0.000000e+00,No           ]  
    ---> setting par type region_size, corresponding to par R
    --->  name = R                 type = region_size               units = cm                val = +1.884609e+17  phys-bounds = [+0.000000e+00,No           ]  
    ---> *** electron distribution parameters ***
    ---> distribution type:  lppl
    ---> r elec. spec. curvature =8.018449e-01
    ---> setting par type curvature, corresponding to par r
    --->  name = r                 type = spectral_curvature        units =                   val = +8.018449e-01  phys-bounds = [-1.000000e+01,+1.000000e+01]  
    ---> s_radio_mm -0.28943913107755925 1.5788782621551185
    ---> s_X 3.9156016040466177
    ---> s_Fermi 1.6729891337010345
    ---> s_UV_X 2.691841699164094
    ---> s_Opt_UV -0.6211800564392291 2.2423601128784583
    ---> s from synch log-log fit -1.0
    ---> s from (s_Fermi + s_UV)/2
    ---> power-law index s, class obj=HSP s chosen is 2.182415
    ---> setting par type LE_spectral_slope, corresponding to par s
    --->  name = s                 type = LE_spectral_slope         units =                   val = +2.182415e+00  phys-bounds = [-1.000000e+01,+1.000000e+01]  
    ---> gamma_3p_Sync= 1.025156e+05, assuming B=5.500000e-02
    ---> gamma_max=1.315786e+06 from nu_max_Sync= 8.544779e+18, using B=5.500000e-02
    ---> setting par type high-energy-cut-off, corresponding to par gmax
    --->  name = gmax              type = high-energy-cut-off       units = Lorentz-factor    val = +1.315786e+06  phys-bounds = [+1.000000e+00,No           ]  
    
    ---> setting par type low-energy-cut-off, corresponding to par gmin
    --->  name = gmin              type = low-energy-cut-off        units = Lorentz-factor    val = +1.350381e+03  phys-bounds = [+1.000000e+00,No           ]  
    
    ---> setting par type turn-over energy, corresponding to par gamma0_log_parab
    ---> using gamma_3p_Sync= 102515.59879590549
    --->  name = gamma0_log_parab  type = turn-over-energy          units = Lorentz-factor    val = +3.169387e+04  phys-bounds = [+1.000000e+00,No           ]  
    
    nu_p_seed_blob 2138672667284302.8
    COMP FACTOR 1.7744345131491372 24949.53232710054
    ---> gamma_3p_SSCc= %e 161476.52800396679
    ---> setting par type turn-over energy, corresponding to par gamma0_log_parab
    ---> using gamma_3p_SSC= 161476.52800396679
    --->  name = gamma0_log_parab  type = turn-over-energy          units = Lorentz-factor    val = +4.992232e+04  phys-bounds = [+1.000000e+00,No           ]  
    
    
    ---> setting par type electron_density, corresponding to par N
    ---> B from nu_p_S=2.216786e-02
    ---> get B from best matching of nu_p_IC
    ---> B=4.012546e-01, out of boundaries 1.000000e-02 1.000000e-01, rejected
         Best B not found, (temporary set to 1.000000e-01)
    ---> setting par type magnetic_field, corresponding to par B
    --->  name = B                 type = magnetic_field            units = G                 val = +1.000000e-01  phys-bounds = [+0.000000e+00,No           ]  
    
    ---> constrain failed, B set to:  name = B                 type = magnetic_field            units = G                 val = +1.000000e-01  phys-bounds = [+0.000000e+00,No           ]  
    
    
    ---> update pars for new B 
    ---> setting par type low-energy-cut-off, corresponding to par gmin
    --->  name = gmin              type = low-energy-cut-off        units = Lorentz-factor    val = +1.001469e+03  phys-bounds = [+1.000000e+00,No           ]  
    
    ---> setting par type low-energy-cut-off, corresponding to par gamma0_log_parab
    ---> using gamma_3p_Sync= 76027.60286939003
    --->  name = gamma0_log_parab  type = turn-over-energy          units = Lorentz-factor    val = +2.350480e+04  phys-bounds = [+1.000000e+00,No           ]  
    ---> gamma_max=9.758134e+05 from nu_max_Sync= 8.544779e+18, using B=1.000000e-01
    ---> setting par type high-energy-cut-off, corresponding to par gmax
    --->  name = gmax              type = high-energy-cut-off       units = Lorentz-factor    val = +9.758134e+05  phys-bounds = [+1.000000e+00,No           ]  
    
    ---> setting par type electron_density, corresponding to par N
    ---> get R from Compoton Dominance (CD)
         Best R=1.592619e+16
    ---> setting par type region_size, corresponding to par R
    --->  name = R                 type = region_size               units = cm                val = +1.592619e+16  phys-bounds = [+0.000000e+00,No           ]  
    
    ---> setting par type electron_density, corresponding to par N
    ---> t_var (days) 0.253519856130165
    
    show pars
    -----------------------------------------------------------------------------------------
    model parameters for jet model:
    electron distribution type = lppl  
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +8.480725e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +1.001469e+03 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +9.758134e+05 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.182415e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +8.018449e-01 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +2.350480e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +1.592619e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
    --------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------
    eval_model
    directory .//lppl_jet_prod/ already existing
    removing existing dir
    the directory .//lppl_jet_prod/ has been created
    =============================================================================================
    


.. code:: ipython3

    constr_Plot=Plot(sed_data,interactive=False)
    jet_model.plot_model(plot_obj=constr_Plot)
    
    constr_Plot.add_data_plot(sed_data,autoscale=True)
    constr_Plot.rescale(y_min=-14,y_max=-8,x_min=8)
    constr_Plot.add_residual_plot(jet_model,sed_data,autoscale=True)


.. parsed-literal::

    directory ./ already existing
    label Sum
    label Sync
    label SSC
    res 0



.. image:: fit_example_files/fit_example_34_1.png


SSC/EC fitting
--------------

.. code:: ipython3

    from jetset.model_manager import FitModel
    from jetset.minimizer import  fit_SED
    
    fit_model=FitModel( jet=jet_model, name='SSC-best-fit',  template=my_shape.host_gal)
    
    fit_model.set('z_cosm','frozen')
    
    fit_model.set('beam_obj','frozen')
    
    fit_model.set('nuFnu_p_host','frozen')
    
    #SEDModel.set('L_host',fit_range=[-10.5,-9.5])
    
    fit_model.show_pars()
        
    best_fit=fit_SED(fit_model,sed_data,10.0**11 ,10**28.0,fitname='SSC-best-fit-lppl',minimizer='leastsqbound')


.. parsed-literal::

    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +8.480725e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +1.001469e+03 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +9.758134e+05 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.182415e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +8.018449e-01 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +2.350480e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +1.592619e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | -1.811391e-02 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    directory ./ already existing
    directory .//SSC-best-fit-lppl/ already existing
    removing existing dir
    the directory .//SSC-best-fit-lppl/ has been created
    directory .//SSC-best-fit-lppl/ already existing
    removing existing dir
    the directory .//SSC-best-fit-lppl/ has been created
    filtering data in fit range = [1.000000e+11,1.000000e+28]
    data length 37
    =============================================================================================
    
    *** start fit process ***
    initial pars: 
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +8.480725e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +1.001469e+03 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +9.758134e+05 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.182415e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +8.018449e-01 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +2.350480e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +1.592619e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +1.000000e-01 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | -1.811391e-02 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    minim function calls=200, res=5.517415, chisq=18.09676res check 5.516280233891562 18.09668593879949
    
    **************************************************************************************************
    Fit report
    
    Model: SSC-best-fit-lppl
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +6.972955e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +8.417369e+02 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +8.567258e+09 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.284001e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +1.487699e+00 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +9.390545e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +2.614888e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +4.844317e-02 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | +1.137794e-01 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=2
    calls=204
    mesg=The relative error between two consecutive iterates is at most 0.000000
    dof=28
    chisq=18.096686, chisq/red=0.646310 null hypothesis sig=0.923688
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     N                | +6.972955e-01 | +1.752621e+05 | +8.480725e-01 | [+0.000000e+00,No           ]
     gmin             | +8.417369e+02 | +1.646526e+08 | +1.001469e+03 | [+1.000000e+00,No           ]
     gmax             | +8.567258e+09 | +2.526656e+15 | +9.758134e+05 | [+1.000000e+00,No           ]
     s                | +2.284001e+00 | +5.676853e-02 | +2.182415e+00 | [-1.000000e+01,+1.000000e+01]
     r                | +1.487699e+00 | +2.602300e-01 | +8.018449e-01 | [-1.000000e+01,+1.000000e+01]
     gamma0_log_parab | +9.390545e+04 | +2.386415e+04 | +2.350480e+04 | [+1.000000e+00,No           ]
     R                | +2.614888e+16 | +5.084047e+15 | +1.592619e+16 | [+0.000000e+00,No           ]
     B                | +4.844317e-02 | +8.533596e-03 | +1.000000e-01 | [+0.000000e+00,No           ]
     beam_obj         | Frozen        | Frozen        | +2.500000e+01 | [+1.000000e+00,No           ]
     z_cosm           | Frozen        | Frozen        | +3.080000e-02 | [+0.000000e+00,No           ]
     nuFnu_p_host     | Frozen        | Frozen        | -1.106522e+01 | [-1.148453e+01,-7.484531e+00]
     nu_scale         | +1.137794e-01 | +4.619238e-02 | -1.811391e-02 | [-5.000000e-01,+5.000000e-01]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    ---> 100000000000.0 1e+28 100
    =============================================================================================
    


.. code:: ipython3

    fit_Plot=Plot(sed_data,interactive=True)
    fit_Plot.add_model_plot(fit_model,label='SSC-best-fit')
    fit_Plot.autoscale()
    fit_Plot.rescale(y_min=-14,y_max=-8,x_min=9.0,x_max=30)
    
    for c in fit_model.components[0].spectral_components:
        fit_Plot.add_model_plot(c.SED,autoscale=False,line_style='--')
    
    for c in fit_model.components:
        fit_Plot.add_model_plot(c.SED,autoscale=False,line_style='--')
        
    fit_Plot.add_data_plot(sed_data,autoscale=False,color='b')
    
    
    fit_Plot.add_residual_plot(jet_model,sed_data)



.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing
    label SSC-best-fit
    label Sum
    label Sync
    label SSC
    label Sum
    label host-galaxy
    res 0



.. image:: fit_example_files/fit_example_37_1.png


.. code:: ipython3

    best_fit=fit_SED(fit_model,sed_data,10.0**11 ,10**28.0,fitname='SSC-best-fit-lppl',minimizer='minuit')


.. parsed-literal::

    directory ./ already existing
    directory .//SSC-best-fit-lppl/ already existing
    removing existing dir
    the directory .//SSC-best-fit-lppl/ has been created
    directory .//SSC-best-fit-lppl/ already existing
    removing existing dir
    the directory .//SSC-best-fit-lppl/ has been created
    filtering data in fit range = [1.000000e+11,1.000000e+28]
    data length 37
    =============================================================================================
    
    *** start fit process ***
    initial pars: 
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +6.972955e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +8.417369e+02 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +8.567258e+09 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.284001e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +1.487699e+00 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +9.390545e+04 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +2.614888e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +4.844317e-02 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | +1.137794e-01 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    dict {'par_0': 0.6972954792580421, 'limit_par_0': (0, None), 'par_1': 841.736934087698, 'limit_par_1': (1, None), 'par_2': 8567257866.993572, 'limit_par_2': (1, None), 'par_3': 2.2840007517058893, 'limit_par_3': (-10, 10), 'par_4': 1.4876990201920481, 'limit_par_4': (-10, 10), 'par_5': 93905.44666422937, 'limit_par_5': (1, None), 'par_6': 2.6148876505417004e+16, 'limit_par_6': (0, None), 'par_7': 0.048443167122354636, 'limit_par_7': (0, None), 'par_8': 0.11377935390949812, 'limit_par_8': (-0.5, 0.5)}
    **************************************************
    *                     MIGRAD                     *
    **************************************************
    
    minim function calls=530, res=6.867547, chisq=17.04726**********************************************************************
    ---------------------------------------------------------------------------------------
    fval = 16.915863106222112 | total call = 535 | ncalls = 524
    edm = 297.9262258169708 (Goal: 1e-05) | up = 1.0
    ---------------------------------------------------------------------------------------
    |          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
    ---------------------------------------------------------------------------------------
    |          False |           True |          False |          False |           True |
    ---------------------------------------------------------------------------------------
    |     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
    ---------------------------------------------------------------------------------------
    |          False |           True |           True |             '' |          False |
    ---------------------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------------------------
    |      | Name  |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
    ----------------------------------------------------------------------------------------------
    |    0 | par_0 =  0.7925  |  0.02105 |          |          |          |  0       |          |
    |    1 | par_1 =  892.3   |  0.0005833 |          |          |          |  0       |          |
    |    2 | par_2 =  8.567E+09 |  0.0005836 |          |          |          |  0       |          |
    |    3 | par_3 =  2.362   |  5.775E-05 |          |          | -10      |  10      |          |
    |    4 | par_4 =  1.673   |  0.0001197 |          |          | -10      |  10      |          |
    |    5 | par_5 =  1.234E+05 |  0.0005833 |          |          |          |  0       |          |
    |    6 | par_6 =  2.615E+16 |  0       |          |          |          |  0       |          |
    |    7 | par_7 =  0.04686 |  4.556E-06 |          |          |          |  0       |          |
    |    8 | par_8 =  0.113   |  0.08101 |          |          | -0.5     |  0.5     |          |
    ----------------------------------------------------------------------------------------------
    
    **********************************************************************
    p (0.7922770044026128, 892.2646829130898, 8567257866.993572, 2.3622270126366516, 1.672913136792486, 123360.304928053, 2.6148876505417004e+16, 0.04685965458374586, 0.11297359969144571)
    res check 5.251571513899901 16.917423368359263
    
    **************************************************************************************************
    Fit report
    
    Model: SSC-best-fit-lppl
    --------------------------------------------------------------------------------------------------------------
    model parameters:
     Name             | Type                     | Units            | value         | phys. boundaries
    --------------------------------------------------------------------------------------------------------------
     N                | electron_density         | cm^-3            | +7.924796e-01 | [+0.000000e+00,No           ]  
     gmin             | low-energy-cut-off       | Lorentz-factor   | +8.922647e+02 | [+1.000000e+00,No           ]  
     gmax             | high-energy-cut-off      | Lorentz-factor   | +8.567258e+09 | [+1.000000e+00,No           ]  
     s                | LE_spectral_slope        |                  | +2.362228e+00 | [-1.000000e+01,+1.000000e+01]  
     r                | spectral_curvature       |                  | +1.672914e+00 | [-1.000000e+01,+1.000000e+01]  
     gamma0_log_parab | turn-over-energy         | Lorentz-factor   | +1.233603e+05 | [+1.000000e+00,No           ]  
     R                | region_size              | cm               | +2.614888e+16 | [+0.000000e+00,No           ]  
     B                | magnetic_field           | G                | +4.685961e-02 | [+0.000000e+00,No           ]  
     beam_obj         | beaming                  |                  | +2.500000e+01 | [+1.000000e+00,No           ]  
     z_cosm           | redshift                 |                  | +3.080000e-02 | [+0.000000e+00,No           ]  
     nuFnu_p_host     | nuFnu-scale              | erg cm^-2 s^-1   | -1.106522e+01 | [-2.000000e+01,+2.000000e+01]  
     nu_scale         | nu-scale                 | Hz               | +1.129707e-01 | [-2.000000e+00,+2.000000e+00]  
    --------------------------------------------------------------------------------------------------------------
    
    converged=True
    calls=536
    mesg=
    dof=28
    chisq=16.917423, chisq/red=0.604194 null hypothesis sig=0.950202
    
    best fit pars
    ---------------------------------------------------------------------------------------------------
    best-fit parameters:
      Name            | best-fit value| best-fit err  | start value   | fit boundaries
    ---------------------------------------------------------------------------------------------------
     N                | +7.924796e-01 | +2.105191e-02 | +6.972955e-01 | [+0.000000e+00,No           ]
     gmin             | +8.922647e+02 | +5.833143e-04 | +8.417369e+02 | [+1.000000e+00,No           ]
     gmax             | +8.567258e+09 | +5.836487e-04 | +8.567258e+09 | [+1.000000e+00,No           ]
     s                | +2.362228e+00 | +5.775267e-05 | +2.284001e+00 | [-1.000000e+01,+1.000000e+01]
     r                | +1.672914e+00 | +1.197320e-04 | +1.487699e+00 | [-1.000000e+01,+1.000000e+01]
     gamma0_log_parab | +1.233603e+05 | +5.833121e-04 | +9.390545e+04 | [+1.000000e+00,No           ]
     R                | +2.614888e+16 | +0.000000e+00 | +2.614888e+16 | [+0.000000e+00,No           ]
     B                | +4.685961e-02 | +4.556030e-06 | +4.844317e-02 | [+0.000000e+00,No           ]
     beam_obj         | Frozen        | Frozen        | +2.500000e+01 | [+1.000000e+00,No           ]
     z_cosm           | Frozen        | Frozen        | +3.080000e-02 | [+0.000000e+00,No           ]
     nuFnu_p_host     | Frozen        | Frozen        | -1.106522e+01 | [-1.148453e+01,-7.484531e+00]
     nu_scale         | +1.129707e-01 | +8.101221e-02 | +1.137794e-01 | [-5.000000e-01,+5.000000e-01]
    ---------------------------------------------------------------------------------------------------
    **************************************************************************************************
    
    ---> 100000000000.0 1e+28 100
    =============================================================================================
    


.. code:: ipython3

    fit_Plot=Plot(sed_data,interactive=True)
    fit_Plot.add_model_plot(fit_model,label='SSC-best-fit')
    fit_Plot.autoscale()
    fit_Plot.rescale(y_min=-14,y_max=-8,x_min=9.0,x_max=30)
    
    for c in fit_model.components[0].spectral_components:
        fit_Plot.add_model_plot(c.SED,autoscale=False,line_style='--')
    
    for c in fit_model.components:
        fit_Plot.add_model_plot(c.SED,autoscale=False,line_style='--')
        
    fit_Plot.add_data_plot(sed_data,autoscale=False,color='b')
    
    
    fit_Plot.add_residual_plot(jet_model,sed_data)
    



.. parsed-literal::

    running PyLab in interactive mode
    directory ./ already existing
    label SSC-best-fit
    label Sum
    label Sync
    label SSC
    label Sum
    label host-galaxy
    res 0



.. image:: fit_example_files/fit_example_39_1.png

