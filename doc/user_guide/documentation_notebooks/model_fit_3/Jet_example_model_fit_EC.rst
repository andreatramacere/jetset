.. _model_fitting_EC:

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')

Model fitting 3: External Compton
=================================

Loading data
------------

see the :ref:`data_format` user guide for further information about loading data and :ref:`jet_physical_guide_EC` for the information regarding the implementation of the external Conpton model

.. code:: ipython3

    from jetset.jet_model import Jet
    from jetset.data_loader import Data,ObsData
    from jetset.test_data_helper import  test_SEDs
    test_SEDs





.. parsed-literal::

    ['/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_3C345.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk421_EBL_DEABS.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_ABS.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_DEABS.ecsv']



.. code:: ipython3

    data=Data.from_file(test_SEDs[0])


.. code:: ipython3

    sed_data=ObsData(data_table=data)

.. code:: ipython3

    %matplotlib inline
    p=sed_data.plot_sed(show_dataset=True)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_8_0.png


we filter out the data set ``-1``

.. code:: ipython3

    sed_data.show_data_sets()
    sed_data.filter_data_set('-1',exclude=True)
    sed_data.show_data_sets()
    p=sed_data.plot_sed()



.. parsed-literal::

    current datasets
    dataset -1
    dataset 0
    dataset 1
    dataset 2
    ---> excluding  data_set/s ['-1']
    filter -1 192
    current datasets
    dataset 0
    dataset 1
    dataset 2
    ---> data sets left after filtering None
    ---> data len after filtering=192
    current datasets
    dataset 0
    dataset 1
    dataset 2



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_10_1.png


.. code:: ipython3

    sed_data.group_data(bin_width=.2)
    sed_data.add_systematics(0.2,[10.**6,10.**29])
    p=sed_data.plot_sed()


.. parsed-literal::

    ================================================================================
    
    ***  binning data  ***
    ---> N bins= 80
    ---> bin_widht= 0.2
    ================================================================================
    



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_11_1.png


.. code:: ipython3

    sed_data.save('3C454_data.pkl')

Phenomenological model constraining
-----------------------------------

see the :ref:`phenom_constr` user guide for further information about phenomenological model constraining

.. code:: ipython3

    from jetset.sed_shaper import  SEDShape
    my_shape=SEDShape(sed_data)
    my_shape.eval_indices(silent=True)
    p=my_shape.plot_indices()
    p.setlim(y_min=1E-15,y_max=1E-6)


.. parsed-literal::

    ================================================================================
    
    *** evaluating spectral indices for data ***
    ================================================================================
    



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_15_1.png


for the synchrotron sed_shaping we include the check for Big Blue Bump
(BBB) component. Moreover, we force the model to use a pure
log-parabolic function and not a log-cubic one in order to get a better
estimation of the BBB component. The fit values of the BBB component
will be used in the ``ObsConstrain`` to guess the accretion disk
luminosity and temperature

.. code:: ipython3

    mm,best_fit=my_shape.sync_fit(check_BBB_template=True,
                                  check_host_gal_template=False,
                                  use_log_par=True,
                                  Ep_start=None,
                                  minimizer='lsb',
                                  silent=True,
                                  fit_range=[9,16])


.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the synchrotron component ***
    ---> first blind fit run,  fit range: [9, 16]
    --> class:  LSP
    
    --> class:  LSP
    
    



.. raw:: html

    <i>Table length=5</i>
    <table id="table140377071359888-703608" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogParabolaEp</td><td>b</td><td>-2.984653e-01</td><td>-2.984653e-01</td><td>5.631694e-02</td><td>--</td><td>-1.527892e-01</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogParabolaEp</td><td>Ep</td><td>1.190850e+01</td><td>1.190850e+01</td><td>2.238841e-01</td><td>--</td><td>1.298338e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogParabolaEp</td><td>Sp</td><td>-1.123366e+01</td><td>-1.123366e+01</td><td>7.306404e-02</td><td>--</td><td>-1.095506e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>BBB</td><td>nuFnu_p_BBB</td><td>-1.155965e+01</td><td>-1.155965e+01</td><td>6.791135e-02</td><td>--</td><td>-1.095506e+01</td><td>-1.295506e+01</td><td>-8.955061e+00</td><td>False</td></tr>
    <tr><td>BBB</td><td>nu_scale</td><td>7.058302e-02</td><td>7.058302e-02</td><td>2.539034e-03</td><td>--</td><td>0.000000e+00</td><td>-5.000000e-01</td><td>5.000000e-01</td><td>False</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377071359888-703608').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377071359888-703608').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> sync       nu_p=+1.190850e+01 (err=+2.238841e-01)  nuFnu_p=-1.123366e+01 (err=+7.306404e-02) curv.=-2.984653e-01 (err=+5.631694e-02)
    ================================================================================
    


.. code:: ipython3

    my_shape.IC_fit(fit_range=[16,26],minimizer='minuit', silent=True)
    p=my_shape.plot_shape_fit()
    p.setlim(y_min=1E-15)


.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the IC component ***
    ---> fit range: [16, 26]
    ---> LogCubic fit
    
    



.. raw:: html

    <i>Table length=4</i>
    <table id="table140377092391792-472312" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-1.128855e-01</td><td>-1.128855e-01</td><td>1.240849e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-1.065003e-02</td><td>-1.065003e-02</td><td>2.393721e-03</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>2.273378e+01</td><td>2.273378e+01</td><td>1.453319e-01</td><td>--</td><td>2.270678e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.043090e+01</td><td>-1.043090e+01</td><td>6.087264e-02</td><td>--</td><td>-1.000000e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377092391792-472312').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377092391792-472312').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> IC         nu_p=+2.273378e+01 (err=+1.453319e-01)  nuFnu_p=-1.043090e+01 (err=+6.087264e-02) curv.=-1.128855e-01 (err=+1.240849e-02)
    ================================================================================
    



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_18_3.png


In this case we use the ``constrain_SSC_EC_model``, and we ask to use a
dusty torus and BLR component external component

read the section :ref:`jet_physical_guide_EC`  for more information regarding the EC model

.. code:: ipython3

    from jetset.obs_constrain import ObsConstrain
    from jetset.model_manager import  FitModel
    from jetset.minimizer import fit_SED
    sed_obspar=ObsConstrain(beaming=25,
                            B_range=[0.1,0.2],
                            distr_e='bkn',
                            t_var_sec=7*86400,
                            nu_cut_IR=1E9,
                            SEDShape=my_shape)
    
    
    prefit_jet=sed_obspar.constrain_SSC_EC_model(electron_distribution_log_values=False,EC_componets_list=['EC_DT','EC_BLR'],R_H=1E18,silent=True)



.. parsed-literal::

    ================================================================================
    
    ***  constrains parameters from observable ***
    



.. raw:: html

    <i>Table length=19</i>
    <table id="table140377057313408-260650" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>2.845488e+17</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.500000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>2.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.071498e+01</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.601124e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.846756e+01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>3.049588e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.357911e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.500000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>DT</td><td>cm</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>BLR</td><td>cm</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>BLR</td><td>cm</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>Disk</td><td>erg / s</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>3.018434e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377057313408-260650').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377057313408-260650').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    
    ================================================================================
    


.. code:: ipython3

    prefit_jet.eval()
    p=prefit_jet.plot_model(sed_data=sed_data)
    prefit_jet.save_model('prefit_jet_EC.pkl')



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_22_0.png


The prefit model should works well for the synchrotron component, but
the EC one is a bit problematic. We can set as starting values a
slightly harder value of ``p``, and a larger value of ``gamma_break``
and ``gmax``. We freeze some parameters, and we also set some
``fit_range`` values. Setting fit_range can speed-up the fit convergence
but should be judged by the user each time according to the physics of
the particular source

.. note::
   With the new implementation of composite model  (`FitModel` class) to set parameters you have to specify the model component, this is different from versions<1.1.2,
   and this holds also for the `freeze` method and for setting  `fit_range` intervals, and for the methods relate to parameters setting in general.
   See the :ref:`composite_models` user guide for further information about the new implementation of `FitModel`, in particular for parameter setting

EC model fit
------------

.. code:: ipython3

    jet=Jet.load_model('prefit_jet_EC.pkl')
    jet.set_gamma_grid_size(100)
    fit_model=FitModel( jet=jet, name='EC-best-fit-lsb')
    fit_model.show_model_components()



.. raw:: html

    <i>Table length=19</i>
    <table id="table140377045339200-197714" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.071498e+01</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.601124e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.846756e+01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>3.049588e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.357911e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.500000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>DT</td><td>cm</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>BLR</td><td>cm</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>BLR</td><td>cm</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>Disk</td><td>erg / s</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>3.018434e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>2.845488e+17</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.500000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>2.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377045339200-197714').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377045339200-197714').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    
    --------------------------------------------------------------------------------
    Composite model description
    --------------------------------------------------------------------------------
    name: EC-best-fit-lsb  
    type: composite_model  
    components models:
     -model name: jet_leptonic model type: jet
    
    --------------------------------------------------------------------------------


.. code:: ipython3

    
    fit_model.freeze('jet_leptonic','z_cosm')
    fit_model.free('jet_leptonic','R_H')
    fit_model.freeze('jet_leptonic','L_Disk')
    fit_model.freeze('jet_leptonic','R_DT')
    fit_model.freeze('jet_leptonic','R_BLR_in')
    fit_model.freeze('jet_leptonic','R_BLR_out')
    
    fit_model.jet_leptonic.parameters.R.fit_range=[1E16,5E18]
    fit_model.jet_leptonic.parameters.gamma_break.fit_range=[300,3000]
    fit_model.jet_leptonic.parameters.gmin.fit_range=[2,100]
    fit_model.jet_leptonic.parameters.gmax.fit_range=[1000,1E6]

.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer_lsb=ModelMinimizer('lsb')
    best_fit_lsb=model_minimizer_lsb.fit(fit_model,sed_data,1E11,1E29,fitname='EC-best-fit-lsb',repeat=3)


.. parsed-literal::

    filtering data in fit range = [1.000000e+11,1.000000e+29]
    data length 21
    ================================================================================
    
    *** start fit process ***
    ----- 
    fit run: 0



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.93395e+01
    
    fit run: 1
    - old chisq=1.93395e+01



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.93395e+01
    
    fit run: 2
    - old chisq=1.93395e+01



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.93386e+01
    
    -------------------------------------------------------------------------
    Fit report
    
    Model: EC-best-fit-lsb



.. raw:: html

    <i>Table length=19</i>
    <table id="table140377067785472-591660" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>4.297117e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>2.822712e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>2.258362e+01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>3.008818e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>1.034138e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.639483e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>4.580689e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>DT</td><td>cm</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.781413e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>2.089862e-04</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>BLR</td><td>cm</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>BLR</td><td>cm</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>Disk</td><td>erg / s</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>3.780580e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>4.119800e+17</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>9.984351e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.277432e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>1.071574e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377067785472-591660').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377067785472-591660').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    
    converged=True
    calls=131
    mesg=



.. parsed-literal::

    'The relative error between two consecutive iterates is at most 0.000000'


.. parsed-literal::

    dof=7
    chisq=19.338558, chisq/red=2.762651 null hypothesis sig=0.007190
    
    best fit pars



.. raw:: html

    <i>Table length=19</i>
    <table id="table140377079496464-101377" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>4.297117e+00</td><td>4.297117e+00</td><td>3.882204e+01</td><td>--</td><td>1.071498e+01</td><td>2.000000e+00</td><td>1.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>2.822712e+04</td><td>2.822712e+04</td><td>6.287440e+04</td><td>--</td><td>1.601124e+04</td><td>1.000000e+03</td><td>1.000000e+06</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>2.258362e+01</td><td>2.258362e+01</td><td>3.495097e+01</td><td>--</td><td>1.846756e+01</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>3.008818e+02</td><td>3.008818e+02</td><td>9.752563e+01</td><td>--</td><td>3.049588e+02</td><td>3.000000e+02</td><td>3.000000e+03</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>1.034138e+00</td><td>1.034138e+00</td><td>8.354114e-01</td><td>--</td><td>2.357911e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>3.639483e+00</td><td>3.639483e+00</td><td>5.168124e-01</td><td>--</td><td>3.500000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>4.580689e+02</td><td>4.580689e+02</td><td>5.222955e+03</td><td>--</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>5.143375e+18</td><td>--</td><td>--</td><td>--</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>1.781413e-01</td><td>1.781413e-01</td><td>3.296905e+00</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>2.089862e-04</td><td>2.089862e-04</td><td>2.564684e+04</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>2.057350e+17</td><td>--</td><td>--</td><td>--</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>4.114700e+17</td><td>--</td><td>--</td><td>--</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>4.232688e+45</td><td>--</td><td>--</td><td>--</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>3.780580e+04</td><td>3.780580e+04</td><td>1.660484e+04</td><td>--</td><td>3.018434e+04</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>4.119800e+17</td><td>4.119800e+17</td><td>1.157586e+18</td><td>--</td><td>2.845488e+17</td><td>1.000000e+16</td><td>5.000000e+18</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>9.984351e+17</td><td>9.984351e+17</td><td>1.675663e+22</td><td>--</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>1.277432e-01</td><td>1.277432e-01</td><td>4.863055e-01</td><td>--</td><td>1.500000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>1.071574e+01</td><td>1.071574e+01</td><td>3.382188e+01</td><td>--</td><td>2.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>5.930000e-01</td><td>--</td><td>--</td><td>--</td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140377079496464-101377').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140377079496464-101377').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    -------------------------------------------------------------------------
    
    ================================================================================
    


.. code:: ipython3

    %matplotlib inline
    fit_model.set_nu_grid(1E6,1E30,200)
    fit_model.eval()
    p2=fit_model.plot_model(sed_data=sed_data)
    p2.setlim(y_min=1E-15,y_max=5E-9,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_29_0.png


.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer_minuit=ModelMinimizer('minuit')
    fit_model.freeze('jet_leptonic','R_H')
    fit_model.jet_leptonic.parameters.gmax.val=1E5
    best_fit_minuit=model_minimizer_minuit.fit(fit_model,sed_data,1E11,1E29,fitname='EC-best-fit-minuit',repeat=2)


.. parsed-literal::

    filtering data in fit range = [1.000000e+11,1.000000e+29]
    data length 21
    ================================================================================
    
    *** start fit process ***
    ----- 
    fit run: 0



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.16376e+01
    
    fit run: 1
    - old chisq=1.16376e+01



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.05577e+01
    
    -------------------------------------------------------------------------
    Fit report
    
    Model: EC-best-fit-minuit



.. raw:: html

    <i>Table length=19</i>
    <table id="table140376602092112-344415" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>4.295675e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.002469e+05</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>2.257972e+01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>3.008410e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>1.036353e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.604263e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>4.431086e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>DT</td><td>cm</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.761708e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>7.299979e-09</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>BLR</td><td>cm</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>BLR</td><td>cm</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>Disk</td><td>erg / s</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>3.853341e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>4.109322e+17</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>9.984351e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.275298e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>1.069429e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140376602092112-344415').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140376602092112-344415').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    
    converged=True
    calls=1240
    mesg=



.. raw:: html

    <table>
        <tr>
            <td colspan="2" style="text-align:left" title="Minimum value of function"> FCN = 10.56 </td>
            <td colspan="3" style="text-align:center" title="No. of function evaluations in last call and total number"> Nfcn = 1240 </td>
        </tr>
        <tr>
            <td colspan="2" style="text-align:left" title="Estimated distance to minimum and goal"> EDM = 9.31e+08 (Goal: 0.0002) </td>
            <td colspan="3" style="text-align:center" title="No. of gradient evaluations in last call and total number">  </td>
        </tr>
        <tr>
            <td style="text-align:center;background-color:#c15ef7;color:black"> INVALID Minimum </td>
            <td style="text-align:center;background-color:#92CCA6;color:black"> Valid Parameters </td>
            <td colspan="3" style="text-align:center;background-color:#92CCA6;color:black"> No Parameters at limit </td>
        </tr>
        <tr>
            <td colspan="2" style="text-align:center;background-color:#c15ef7;color:black"> ABOVE EDM threshold (goal x 10) </td>
            <td colspan="3" style="text-align:center;background-color:#92CCA6;color:black"> Below call limit </td>
        </tr>
        <tr>
            <td style="text-align:center;background-color:#92CCA6;color:black"> Covariance </td>
            <td style="text-align:center;background-color:#92CCA6;color:black"> Hesse ok </td>
            <td style="text-align:center;background-color:#FFF79A;color:black" title="Is covariance matrix accurate?"> APPROXIMATE </td>
            <td style="text-align:center;background-color:#c15ef7;color:black" title="Is covariance matrix positive definite?"> NOT pos. def. </td>
            <td style="text-align:center;background-color:#c15ef7;color:black" title="Was positive definiteness enforced by Minuit?"> FORCED </td>
        </tr>
    </table><table>
        <tr>
            <td></td>
            <th title="Variable name"> Name </th>
            <th title="Value of parameter"> Value </th>
            <th title="Hesse error"> Hesse Error </th>
            <th title="Minos lower error"> Minos Error- </th>
            <th title="Minos upper error"> Minos Error+ </th>
            <th title="Lower limit of the parameter"> Limit- </th>
            <th title="Upper limit of the parameter"> Limit+ </th>
            <th title="Is the parameter fixed in the fit"> Fixed </th>
        </tr>
        <tr>
            <th> 0 </th>
            <td> par_0 </td>
            <td> 4.2956749 </td>
            <td> 0.0000006 </td>
            <td>  </td>
            <td>  </td>
            <td> 2 </td>
            <td> 100 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 1 </th>
            <td> par_1 </td>
            <td> 100.246922e3 </td>
            <td> 0.000015e3 </td>
            <td>  </td>
            <td>  </td>
            <td> 1E+03 </td>
            <td> 1E+06 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 2 </th>
            <td> par_2 </td>
            <td> 22.5797225 </td>
            <td> 0.0000005 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td>  </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 3 </th>
            <td> par_3 </td>
            <td> 300.841009 </td>
            <td> 0.000011 </td>
            <td>  </td>
            <td>  </td>
            <td> 300 </td>
            <td> 3E+03 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 4 </th>
            <td> par_4 </td>
            <td> 1 </td>
            <td> 6 </td>
            <td>  </td>
            <td>  </td>
            <td> -10 </td>
            <td> 10 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 5 </th>
            <td> par_5 </td>
            <td> 3.604263 </td>
            <td> 0.000007 </td>
            <td>  </td>
            <td>  </td>
            <td> -10 </td>
            <td> 10 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 6 </th>
            <td> par_6 </td>
            <td> 443.10861749 </td>
            <td> 0.00000004 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td>  </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 7 </th>
            <td> par_7 </td>
            <td> 176.170830e-3 </td>
            <td> 0.000015e-3 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td> 1 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 8 </th>
            <td> par_8 </td>
            <td> 7.3000e-9 </td>
            <td> 0.0034e-9 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td> 1 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 9 </th>
            <td> par_9 </td>
            <td> 38.53341003833e3 </td>
            <td> 0.00000000004e3 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td>  </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 10 </th>
            <td> par_10 </td>
            <td> 410.93223e15 </td>
            <td> 0.00027e15 </td>
            <td>  </td>
            <td>  </td>
            <td> 1E+16 </td>
            <td> 5E+18 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 11 </th>
            <td> par_11 </td>
            <td> 127.529767e-3 </td>
            <td> 0.000023e-3 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td>  </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 12 </th>
            <td> par_12 </td>
            <td> 10.69428574 </td>
            <td> 0.00000004 </td>
            <td>  </td>
            <td>  </td>
            <td> 0.0001 </td>
            <td>  </td>
            <td>  </td>
        </tr>
    </table><table>
        <tr>
            <td></td>
            <th> par_0 </th>
            <th> par_1 </th>
            <th> par_2 </th>
            <th> par_3 </th>
            <th> par_4 </th>
            <th> par_5 </th>
            <th> par_6 </th>
            <th> par_7 </th>
            <th> par_8 </th>
            <th> par_9 </th>
            <th> par_10 </th>
            <th> par_11 </th>
            <th> par_12 </th>
        </tr>
        <tr>
            <th> par_0 </th>
            <td> 4.1e-13 </td>
            <td style="background-color:rgb(250,218,218);color:black"> 2.04e-09 <strong>(0.215)</strong> </td>
            <td style="background-color:rgb(250,196,196);color:black"> 1.21e-13 <strong>(0.363)</strong> </td>
            <td style="background-color:rgb(250,196,196);color:black"> 2.62e-12 <strong>(0.359)</strong> </td>
            <td style="background-color:rgb(250,195,195);color:black"> 1.42e-06 <strong>(0.364)</strong> </td>
            <td style="background-color:rgb(250,195,195);color:black"> 1.51e-12 <strong>(0.363)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.12e-22 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 4.92e-19 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.75e-27 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.69e-24 </td>
            <td style="background-color:rgb(250,197,197);color:black"> 6.07e+04 <strong>(0.356)</strong> </td>
            <td style="background-color:rgb(222,222,250);color:black"> -3.23e-15 <strong>(-0.218)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.26e-18 </td>
        </tr>
        <tr>
            <th> par_1 </th>
            <td style="background-color:rgb(250,218,218);color:black"> 2.04e-09 <strong>(0.215)</strong> </td>
            <td> 0.000221 </td>
            <td style="background-color:rgb(250,162,162);color:black"> 4.54e-09 <strong>(0.587)</strong> </td>
            <td style="background-color:rgb(250,163,163);color:black"> 9.84e-08 <strong>(0.581)</strong> </td>
            <td style="background-color:rgb(250,162,162);color:black"> 0.0533 <strong>(0.589)</strong> </td>
            <td style="background-color:rgb(250,162,162);color:black"> 5.69e-08 <strong>(0.588)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.93e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.85e-14 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.16e-22 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.37e-20 </td>
            <td style="background-color:rgb(250,164,164);color:black"> 2.28e+09 <strong>(0.577)</strong> </td>
            <td style="background-color:rgb(204,204,250);color:black"> -1.22e-10 <strong>(-0.353)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -8.49e-14 </td>
        </tr>
        <tr>
            <th> par_2 </th>
            <td style="background-color:rgb(250,196,196);color:black"> 1.21e-13 <strong>(0.363)</strong> </td>
            <td style="background-color:rgb(250,162,162);color:black"> 4.54e-09 <strong>(0.587)</strong> </td>
            <td> 2.7e-13 </td>
            <td style="background-color:rgb(250,103,103);color:black"> 5.81e-12 <strong>(0.983)</strong> </td>
            <td style="background-color:rgb(250,100,100);color:black"> 3.15e-06 <strong>(0.997)</strong> </td>
            <td style="background-color:rgb(250,101,101);color:black"> 3.36e-12 <strong>(0.995)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.14e-21 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.09e-18 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.28e-26 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -3.76e-24 </td>
            <td style="background-color:rgb(250,104,104);color:black"> 1.35e+05 <strong>(0.976)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -7.19e-15 <strong>(-0.597)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.02e-18 </td>
        </tr>
        <tr>
            <th> par_3 </th>
            <td style="background-color:rgb(250,196,196);color:black"> 2.62e-12 <strong>(0.359)</strong> </td>
            <td style="background-color:rgb(250,163,163);color:black"> 9.84e-08 <strong>(0.581)</strong> </td>
            <td style="background-color:rgb(250,103,103);color:black"> 5.81e-12 <strong>(0.983)</strong> </td>
            <td> 1.3e-10 </td>
            <td style="background-color:rgb(250,102,102);color:black"> 6.83e-05 <strong>(0.986)</strong> </td>
            <td style="background-color:rgb(250,102,102);color:black"> 7.29e-11 <strong>(0.984)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.47e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.37e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.77e-25 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -8.15e-23 </td>
            <td style="background-color:rgb(250,105,105);color:black"> 2.92e+06 <strong>(0.965)</strong> </td>
            <td style="background-color:rgb(173,173,250);color:black"> -1.56e-13 <strong>(-0.590)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.09e-16 </td>
        </tr>
        <tr>
            <th> par_4 </th>
            <td style="background-color:rgb(250,195,195);color:black"> 1.42e-06 <strong>(0.364)</strong> </td>
            <td style="background-color:rgb(250,162,162);color:black"> 0.0533 <strong>(0.589)</strong> </td>
            <td style="background-color:rgb(250,100,100);color:black"> 3.15e-06 <strong>(0.997)</strong> </td>
            <td style="background-color:rgb(250,102,102);color:black"> 6.83e-05 <strong>(0.986)</strong> </td>
            <td> 37 </td>
            <td style="background-color:rgb(250,100,100);color:black"> 3.95e-05 <strong>(0.998)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.34e-14 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.29e-11 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.5e-19 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.42e-17 </td>
            <td style="background-color:rgb(250,103,103);color:black"> 1.58e+12 <strong>(0.979)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -8.45e-08 <strong>(-0.599)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.9e-11 </td>
        </tr>
        <tr>
            <th> par_5 </th>
            <td style="background-color:rgb(250,195,195);color:black"> 1.51e-12 <strong>(0.363)</strong> </td>
            <td style="background-color:rgb(250,162,162);color:black"> 5.69e-08 <strong>(0.588)</strong> </td>
            <td style="background-color:rgb(250,101,101);color:black"> 3.36e-12 <strong>(0.995)</strong> </td>
            <td style="background-color:rgb(250,102,102);color:black"> 7.29e-11 <strong>(0.984)</strong> </td>
            <td style="background-color:rgb(250,100,100);color:black"> 3.95e-05 <strong>(0.998)</strong> </td>
            <td> 4.23e-11 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.43e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.37e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.6e-25 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.72e-23 </td>
            <td style="background-color:rgb(250,103,103);color:black"> 1.69e+06 <strong>(0.977)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -9.01e-14 <strong>(-0.597)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.29e-17 </td>
        </tr>
        <tr>
            <th> par_6 </th>
            <td style="background-color:rgb(250,250,250);color:black"> -5.12e-22 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.93e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.14e-21 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.47e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.34e-14 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.43e-20 </td>
            <td> 1.62e-15 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.64e-27 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 5.42e-35 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.6e-32 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -0.000572 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 3.05e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.13e-26 </td>
        </tr>
        <tr>
            <th> par_7 </th>
            <td style="background-color:rgb(250,250,250);color:black"> 4.92e-19 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.85e-14 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.09e-18 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.37e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.29e-11 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.37e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.64e-27 </td>
            <td> 2.35e-16 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.21e-32 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.53e-29 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 0.55 <strong>(0.000)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.93e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.05e-23 </td>
        </tr>
        <tr>
            <th> par_8 </th>
            <td style="background-color:rgb(250,250,250);color:black"> -5.75e-27 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.16e-22 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.28e-26 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.77e-25 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.5e-19 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.6e-25 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 5.42e-35 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.21e-32 </td>
            <td> 1.18e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.79e-37 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.42e-09 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 3.42e-28 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.39e-31 </td>
        </tr>
        <tr>
            <th> par_9 </th>
            <td style="background-color:rgb(250,250,250);color:black"> -1.69e-24 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.37e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -3.76e-24 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -8.15e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.42e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -4.72e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.6e-32 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.53e-29 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.79e-37 </td>
            <td> 1.62e-15 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.89e-06 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.01e-25 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 7.04e-29 </td>
        </tr>
        <tr>
            <th> par_10 </th>
            <td style="background-color:rgb(250,197,197);color:black"> 6.07e+04 <strong>(0.356)</strong> </td>
            <td style="background-color:rgb(250,164,164);color:black"> 2.28e+09 <strong>(0.577)</strong> </td>
            <td style="background-color:rgb(250,104,104);color:black"> 1.35e+05 <strong>(0.976)</strong> </td>
            <td style="background-color:rgb(250,105,105);color:black"> 2.92e+06 <strong>(0.965)</strong> </td>
            <td style="background-color:rgb(250,103,103);color:black"> 1.58e+12 <strong>(0.979)</strong> </td>
            <td style="background-color:rgb(250,103,103);color:black"> 1.69e+06 <strong>(0.977)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -0.000572 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 0.55 <strong>(0.000)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.42e-09 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.89e-06 </td>
            <td> 7.08e+22 </td>
            <td style="background-color:rgb(174,174,250);color:black"> -3.61e+03 <strong>(-0.586)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.52 <strong>(-0.000)</strong> </td>
        </tr>
        <tr>
            <th> par_11 </th>
            <td style="background-color:rgb(222,222,250);color:black"> -3.23e-15 <strong>(-0.218)</strong> </td>
            <td style="background-color:rgb(204,204,250);color:black"> -1.22e-10 <strong>(-0.353)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -7.19e-15 <strong>(-0.597)</strong> </td>
            <td style="background-color:rgb(173,173,250);color:black"> -1.56e-13 <strong>(-0.590)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -8.45e-08 <strong>(-0.599)</strong> </td>
            <td style="background-color:rgb(172,172,250);color:black"> -9.01e-14 <strong>(-0.597)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> 3.05e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.93e-20 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 3.42e-28 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.01e-25 </td>
            <td style="background-color:rgb(174,174,250);color:black"> -3.61e+03 <strong>(-0.586)</strong> </td>
            <td> 5.38e-16 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.35e-19 </td>
        </tr>
        <tr>
            <th> par_12 </th>
            <td style="background-color:rgb(250,250,250);color:black"> -2.26e-18 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -8.49e-14 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.02e-18 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -1.09e-16 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -5.9e-11 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -6.29e-17 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.13e-26 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.05e-23 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 2.39e-31 </td>
            <td style="background-color:rgb(250,250,250);color:black"> 7.04e-29 </td>
            <td style="background-color:rgb(250,250,250);color:black"> -2.52 <strong>(-0.000)</strong> </td>
            <td style="background-color:rgb(250,250,250);color:black"> 1.35e-19 </td>
            <td> 1.61e-15 </td>
        </tr>
    </table>


.. parsed-literal::

    dof=8
    chisq=10.557743, chisq/red=1.319718 null hypothesis sig=0.228038
    
    best fit pars



.. raw:: html

    <i>Table length=19</i>
    <table id="table140376602206608-307236" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>4.295675e+00</td><td>4.295675e+00</td><td>6.400531e-07</td><td>--</td><td>4.297117e+00</td><td>2.000000e+00</td><td>1.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>1.002469e+05</td><td>1.002469e+05</td><td>1.487448e-02</td><td>--</td><td>1.000000e+05</td><td>1.000000e+03</td><td>1.000000e+06</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>2.257972e+01</td><td>2.257972e+01</td><td>5.195583e-07</td><td>--</td><td>2.258362e+01</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>3.008410e+02</td><td>3.008410e+02</td><td>1.138513e-05</td><td>--</td><td>3.008818e+02</td><td>3.000000e+02</td><td>3.000000e+03</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>1.036353e+00</td><td>1.036353e+00</td><td>5.713136e+00</td><td>--</td><td>1.034138e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>3.604263e+00</td><td>3.604263e+00</td><td>6.503548e-06</td><td>--</td><td>3.639483e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>4.431086e+02</td><td>4.431086e+02</td><td>4.021774e-08</td><td>--</td><td>4.580689e+02</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_DT</td><td>5.143375e+18</td><td>--</td><td>--</td><td>--</td><td>5.143375e+18</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>1.761708e-01</td><td>1.761708e-01</td><td>1.532161e-08</td><td>--</td><td>1.781413e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>7.299979e-09</td><td>7.299979e-09</td><td>3.436224e-12</td><td>--</td><td>2.089862e-04</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_in</td><td>2.057350e+17</td><td>--</td><td>--</td><td>--</td><td>2.057350e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_BLR_out</td><td>4.114700e+17</td><td>--</td><td>--</td><td>--</td><td>4.114700e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk</td><td>4.232688e+45</td><td>--</td><td>--</td><td>--</td><td>4.232688e+45</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>3.853341e+04</td><td>3.853341e+04</td><td>4.021422e-08</td><td>--</td><td>3.780580e+04</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>4.109322e+17</td><td>4.109322e+17</td><td>2.660756e+11</td><td>--</td><td>4.119800e+17</td><td>1.000000e+16</td><td>5.000000e+18</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>9.984351e+17</td><td>--</td><td>--</td><td>--</td><td>9.984351e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>1.275298e-01</td><td>1.275298e-01</td><td>2.319318e-08</td><td>--</td><td>1.277432e-01</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>1.069429e+01</td><td>1.069429e+01</td><td>4.007054e-08</td><td>--</td><td>1.071574e+01</td><td>1.000000e-04</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>5.930000e-01</td><td>--</td><td>--</td><td>--</td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140376602206608-307236').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140376602206608-307236').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    -------------------------------------------------------------------------
    
    ================================================================================
    


.. code:: ipython3

    %matplotlib inline
    fit_model.set_nu_grid(1E6,1E30,500)
    fit_model.eval()
    p2=fit_model.plot_model(sed_data=sed_data)
    p2.setlim(y_min=1E-15,y_max=5E-9,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_31_0.png


.. code:: ipython3

    jet.energetic_report()



.. raw:: html

    <i>Table length=37</i>
    <table id="table140376582213984-793546" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
    <tr><td>BulkLorentzFactor</td><td></td><td></td><td>1.069429e+01</td></tr>
    <tr><td>U_e</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.829260e-03</td></tr>
    <tr><td>U_p_cold</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>3.394356e-03</td></tr>
    <tr><td>U_B</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>6.471177e-04</td></tr>
    <tr><td>U_p</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_p_target</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Synch</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.927297e-04</td></tr>
    <tr><td>U_Synch_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>2.520901e+00</td></tr>
    <tr><td>U_Disk</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.331209e-04</td></tr>
    <tr><td>U_BLR</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>3.467122e-11</td></tr>
    <tr><td>U_DT</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.718776e-02</td></tr>
    <tr><td>U_CMB</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Disk_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>1.112494e-02</td></tr>
    <tr><td>U_BLR_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>8.489284e-11</td></tr>
    <tr><td>U_DT_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>7.530747e-05</td></tr>
    <tr><td>U_CMB_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>L_Sync_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>4.051976e+42</td></tr>
    <tr><td>L_SSC_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>9.051115e+41</td></tr>
    <tr><td>L_EC_Disk_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_BLR_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>1.863765e+35</td></tr>
    <tr><td>L_EC_DT_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>1.030114e+44</td></tr>
    <tr><td>L_EC_CMB_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_pp_gamma_rf</td><td>Lum. blob rest. frme.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_Sync</td><td>jet Lum.</td><td>erg / s</td><td>1.158539e+44</td></tr>
    <tr><td>jet_L_SSC</td><td>jet Lum.</td><td>erg / s</td><td>2.587889e+43</td></tr>
    <tr><td>jet_L_EC_Disk</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_BLR</td><td>jet Lum.</td><td>erg / s</td><td>5.328865e+36</td></tr>
    <tr><td>jet_L_EC_DT</td><td>jet Lum.</td><td>erg / s</td><td>2.945295e+45</td></tr>
    <tr><td>jet_L_EC_CMB</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_pp_gamma</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_rad</td><td>jet Lum.</td><td>erg / s</td><td>3.087028e+45</td></tr>
    <tr><td>jet_L_kin</td><td>jet Lum.</td><td>erg / s</td><td>9.459731e+45</td></tr>
    <tr><td>jet_L_tot</td><td>jet Lum.</td><td>erg / s</td><td>1.371866e+46</td></tr>
    <tr><td>jet_L_e</td><td>jet Lum.</td><td>erg / s</td><td>3.312706e+45</td></tr>
    <tr><td>jet_L_B</td><td>jet Lum.</td><td>erg / s</td><td>1.171901e+45</td></tr>
    <tr><td>jet_L_p_cold</td><td>jet Lum.</td><td>erg / s</td><td>6.147025e+45</td></tr>
    <tr><td>jet_L_p</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    
    var astropy_sort_num = function(a, b) {
        var a_num = parseFloat(a);
        var b_num = parseFloat(b);
    
        if (isNaN(a_num) && isNaN(b_num))
            return ((a < b) ? -1 : ((a > b) ? 1 : 0));
        else if (!isNaN(a_num) && !isNaN(b_num))
            return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
        else
            return isNaN(a_num) ? -1 : 1;
    }
    
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table140376582213984-793546').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140376582213984-793546').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [3], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    best_fit_minuit.save_report('EC-best-fit-minuit.pkl')
    model_minimizer_minuit.save_model('EC_model_minimizer_minuit.pkl')
    fit_model.save_model('EC_fit_model_minuit.pkl')

MCMC
----

.. code:: ipython3

    from jetset.mcmc import McmcSampler
    from jetset.minimizer import ModelMinimizer
    model_minimizer_minuit = ModelMinimizer.load_model('EC_model_minimizer_minuit.pkl')


We use a flat prior centered on the best fit value. Setting
``bound=5.0`` and ``bound_rel=True`` means that:

1) the prior interval will be defined as [best_fit_val - delta_m ,
   best_fit_val + delta_p]

2) with delta_p=delta_m=best_fit_val*bound

If we set ``bound_rel=False`` then delta_p = delta_m =
best_fit_err*bound

It is possible to define asymmetric boundaries e.g. ``bound=[2.0,5.0]``
meaning that

1) for ``bound_rel=True``

   delta_p = best_fit_val*bound[1]

   delta_m =b est_fit_val*bound[0]

2) for ``bound_rel=False``

   delta_p = best_fit_err*bound[1]

   delta_m = best_fit_err*bound[0]

In the next release a more flexible prior interface will be added,
including different type of priors

Given the large parameter space, we select a sub sample of parameters
using the ``use_labels_dict``. If we do not pass the ‘use_labels_dict’
the full set of free parameters will be used

.. code:: ipython3

    mcmc=McmcSampler(model_minimizer_minuit)
    
    labels=['N','B','beam_obj','p_1','gamma_break']
    model_name='jet_leptonic'
    use_labels_dict={model_name:labels}
    
    mcmc.run_sampler(nwalkers=128,burnin=10,steps=50,bound=5.0,bound_rel=True,threads=None,walker_start_bound=0.005,use_labels_dict=use_labels_dict)


.. parsed-literal::

    mcmc run starting
    



.. parsed-literal::

      0%|          | 0/50 [00:00<?, ?it/s]


.. parsed-literal::

    mcmc run done, with 1 threads took 267.81 seconds


.. code:: ipython3

    print(mcmc.acceptance_fraction)


.. parsed-literal::

    0.5153125000000001


.. code:: ipython3

    mcmc.model.set_nu_grid(1E6,1E30,200)
    
    p=mcmc.plot_model(sed_data=sed_data,fit_range=[1E11, 1E27],size=50)
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_39_0.png


.. code:: ipython3

    f=mcmc.corner_plot()



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_40_0.png


.. code:: ipython3

    f=mcmc.plot_chain('p_1',log_plot=False)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_41_0.png


Save and reuse MCMC
-------------------

.. code:: ipython3

    mcmc.save('mcmc_sampler.pkl')

.. code:: ipython3

    from jetset.mcmc import McmcSampler
    from jetset.data_loader import ObsData
    from jetset.plot_sedfit import PlotSED
    from jetset.test_data_helper import  test_SEDs
    
    sed_data=ObsData.load('3C454_data.pkl')
    
    ms=McmcSampler.load('mcmc_sampler.pkl')

.. code:: ipython3

    ms.model.set_nu_grid(1E6,1E30,200)
    
    p=ms.plot_model(sed_data=sed_data,fit_range=[1E11, 1E27],size=50)
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_45_0.png


.. code:: ipython3

    f=ms.plot_par('beam_obj',log_plot=False)




.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_46_0.png


.. code:: ipython3

    f=ms.corner_plot()



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_47_0.png


