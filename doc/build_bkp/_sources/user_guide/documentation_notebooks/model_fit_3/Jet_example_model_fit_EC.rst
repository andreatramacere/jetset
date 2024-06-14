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

    ['/Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/test_data/SEDs_data/SED_3C345.ecsv',
     '/Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk421_EBL_DEABS.ecsv',
     '/Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_ABS.ecsv',
     '/Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_DEABS.ecsv']



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
    sed_data.filter_data_set('2',exclude=True)
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
    ---> excluding  data_set/s ['2']
    filter 2 191
    current datasets
    dataset 0
    dataset 1
    ---> data sets left after filtering None
    ---> data len after filtering=191
    current datasets
    dataset 0
    dataset 1



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_10_1.png


.. code:: ipython3

    sed_data.group_data(bin_width=.15)
    sed_data.add_systematics(0.1,[10.**6,10.**29])
    #sed_data.add_systematics(0.05,[10.**19,10.**30])
    
    p=sed_data.plot_sed()


.. parsed-literal::

    ================================================================================
    
    ***  binning data  ***
    ---> N bins= 100
    ---> bin_widht= 0.15
    msk [False  True False  True  True  True  True False  True  True  True  True
      True False  True  True  True False  True False False False False False
     False False False False False False False False False False False False
      True  True  True  True False False False False False False False False
     False False False  True  True  True  True  True  True  True  True False
     False False False False False False False False False False False False
     False False False False False False False False False False False False
     False False False False  True False False  True False False False  True
     False False  True False]
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
    p.setlim(y_min=1E-15,y_max=1E-9)


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
    <table id="table5573765936-845517" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogParabolaEp</td><td>b</td><td>-3.083929e-01</td><td>-3.083929e-01</td><td>2.309292e-02</td><td>--</td><td>-1.560623e-01</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogParabolaEp</td><td>Ep</td><td>1.168757e+01</td><td>1.168757e+01</td><td>9.352232e-02</td><td>--</td><td>1.281229e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogParabolaEp</td><td>Sp</td><td>-1.122282e+01</td><td>-1.122282e+01</td><td>3.377790e-02</td><td>--</td><td>-1.089599e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>BBB</td><td>nuFnu_p_BBB</td><td>-1.155675e+01</td><td>-1.155675e+01</td><td>1.991849e-02</td><td>--</td><td>-1.089599e+01</td><td>-1.289599e+01</td><td>-8.895992e+00</td><td>False</td></tr>
    <tr><td>BBB</td><td>nu_scale</td><td>9.430361e-03</td><td>9.430361e-03</td><td>4.189487e-04</td><td>--</td><td>0.000000e+00</td><td>-5.000000e-01</td><td>5.000000e-01</td><td>False</td></tr>
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
        console.log("$('#table5573765936-845517').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5573765936-845517').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> sync       nu_p=+1.168757e+01 (err=+9.352232e-02)  nuFnu_p=-1.122282e+01 (err=+3.377790e-02) curv.=-3.083929e-01 (err=+2.309292e-02)
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
    ====> simplex
    ====> migrad
    ====> simplex
    ====> migrad
    ====> simplex
    ====> migrad
    
    



.. raw:: html

    <i>Table length=4</i>
    <table id="table5577467664-595076" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-1.239523e-01</td><td>-1.239523e-01</td><td>1.307983e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-1.215669e-02</td><td>-1.215669e-02</td><td>2.221392e-03</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>2.238469e+01</td><td>2.238469e+01</td><td>1.156381e-01</td><td>--</td><td>2.235747e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.038704e+01</td><td>-1.038704e+01</td><td>4.359969e-02</td><td>--</td><td>-1.000000e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table5577467664-595076').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5577467664-595076').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> IC         nu_p=+2.238469e+01 (err=+1.156381e-01)  nuFnu_p=-1.038704e+01 (err=+4.359969e-02) curv.=-1.239523e-01 (err=+1.307983e-02)
    ================================================================================
    



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_18_3.png


In this case we use the ``constrain_SSC_EC_model``, and we ask to use a
dusty torus and BLR component external component

read the section :ref:`jet_physical_guide_EC`  for more information regarding the EC model

.. code:: ipython3

    from jetset.obs_constrain import ObsConstrain
    from jetset.minimizer import fit_SED
    sed_obspar=ObsConstrain(B_range=[0.1,0.2],
                            distr_e='bkn',
                            t_var_sec=15*86400,
                            nu_cut_IR=1E9,
                            theta=2,
                            bulk_factor=20,
                            SEDShape=my_shape)
    
    
    prefit_jet=sed_obspar.constrain_SSC_EC_model(electron_distribution_log_values=False,EC_components_list=['EC_DT','EC_BLR'],R_H=2E18,silent=True,)



.. parsed-literal::

    ================================================================================
    
    ***  constrains parameters from observable ***
    
    ===> setting C threads to 12
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5



.. raw:: html

    <i>Table length=21</i>
    <table id="table5554635872-657504" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>7.607512e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>2.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.500000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>2.000000e+00</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>2.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.033091e+01</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.311585e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>8.405479e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>2.279937e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.309926e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.500000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>1.305539e+19</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.958309e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.154140e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>4.261082e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table5554635872-657504').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5554635872-657504').dataTable({
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



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_22_0.png


.. code:: ipython3

    prefit_jet.make_conical_jet(theta_open=5)


.. parsed-literal::

    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H
    setting R_H to 8.695425919436543e+17


.. code:: ipython3

    prefit_jet.set_EC_dependencies()


.. parsed-literal::

    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5


.. code:: ipython3

    prefit_jet.set_external_field_transf('disk')

.. code:: ipython3

    prefit_jet.eval()
    p=prefit_jet.plot_model(sed_data=sed_data)
    prefit_jet.save_model('prefit_jet_EC.pkl')




.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_26_0.png


The prefit model should works well for the synchrotron component, but
the EC one is a bit problematic. We can set as starting values a
slightly harder value of ``p``, and a larger value of ``gamma_break``
and ``gmax``. We freeze some parameters, and we also set some
``fit_range`` values. Setting fit_range can speed-up the fit convergence
but should be judged by the user each time according to the physics of
the particular source

EC model fit
------------

.. note::
    Please, read the introduction and the caveat :ref:`for the frequentist model fitting <frequentist_model_fitting>` to understand the frequentist fitting workflow
    see the :ref:`composite_models` user guide for further information about the implementation of :class:`.FitModel`, in particular for parameter setting

.. code:: ipython3

    from jetset.data_loader import ObsData
    sed_data=ObsData.load('3C454_data.pkl')
    from jetset.jet_model import Jet

.. code:: ipython3

    from jetset.model_manager import  FitModel
    jet=Jet.load_model('prefit_jet_EC.pkl')
    jet.set_gamma_grid_size(100)
    fit_model=FitModel( jet=jet, name='EC-best-fit-lsb')
    fit_model.show_model_components()


.. parsed-literal::

    ===> setting C threads to 12
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H
    
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
    fit_model.freeze('jet_leptonic','theta')
    
    fit_model.free('jet_leptonic','R_H')
    fit_model.freeze('jet_leptonic','L_Disk')
    fit_model.freeze('jet_leptonic','tau_DT')
    fit_model.freeze('jet_leptonic','tau_BLR')
    
    fit_model.jet_leptonic.parameters.R_H.fit_range=[8E17,5E19]
    fit_model.jet_leptonic.parameters.T_Disk.fit_range=[1E4,1E5]
    fit_model.jet_leptonic.parameters.T_DT.fit_range=[100,1000]
    fit_model.jet_leptonic.parameters.gamma_break.fit_range=[100,500]
    fit_model.jet_leptonic.parameters.gmin.fit_range=[2,100]
    fit_model.jet_leptonic.parameters.gmax.fit_range=[1E4,1E5]
    fit_model.jet_leptonic.parameters.B.fit_range=[1E-2,1]
    fit_model.jet_leptonic.parameters.p.fit_range=[1,2.5]
    fit_model.jet_leptonic.parameters.p_1.fit_range=[3,4]
    fit_model.jet_leptonic.parameters.theta_open.fit_range=[4,6]
    fit_model.jet_leptonic.parameters.BulkFactor.fit_range=[10,30]

.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer_lsb=ModelMinimizer('lsb')
    best_fit_lsb=model_minimizer_lsb.fit(fit_model,sed_data,3E10,1E29,fitname='EC-best-fit-lsb',repeat=1)


.. parsed-literal::

    filtering data in fit range = [3.000000e+10,1.000000e+29]
    data length 24
    ================================================================================
    
    *** start fit process ***
    ----- 



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    - best chisq=1.71336e+02
    
    -------------------------------------------------------------------------
    Fit report
    
    Model: EC-best-fit-lsb



.. raw:: html

    <i>Table length=22</i>
    <table id="table5584434080-65376" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>6.675472e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.581028e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>5.485583e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.894755e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>1.588509e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.466529e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>4.502140e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>1.305539e+19</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.958309e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.154140e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>4.261082e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>2.991580e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>9.134742e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>8.694991e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>9.141046e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>2.000000e+00</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.842049e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.997353e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table5584434080-65376').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5584434080-65376').dataTable({
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
    calls=633
    mesg=



.. parsed-literal::

    '`ftol` termination condition is satisfied.'


.. parsed-literal::

    dof=12
    chisq=171.335707, chisq/red=14.277976 null hypothesis sig=0.000000
    
    best fit pars



.. raw:: html

    <i>Table length=22</i>
    <table id="table5594993216-572158" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>6.675472e+00</td><td>6.675472e+00</td><td>2.084490e+01</td><td>--</td><td>1.033091e+01</td><td>2.000000e+00</td><td>1.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>1.581028e+04</td><td>1.581028e+04</td><td>1.907208e+04</td><td>--</td><td>1.311585e+04</td><td>1.000000e+04</td><td>1.000000e+05</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>5.485583e+02</td><td>5.485583e+02</td><td>4.158130e+03</td><td>--</td><td>8.405479e+02</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>1.894755e+02</td><td>1.894755e+02</td><td>3.075634e+02</td><td>--</td><td>2.279937e+02</td><td>1.000000e+02</td><td>5.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>1.588509e+00</td><td>1.588509e+00</td><td>2.908949e+00</td><td>--</td><td>2.309926e+00</td><td>1.000000e+00</td><td>2.500000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>3.466529e+00</td><td>3.466529e+00</td><td>3.550985e-01</td><td>--</td><td>3.500000e+00</td><td>3.000000e+00</td><td>4.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>4.502140e+02</td><td>4.502140e+02</td><td>3.219246e+03</td><td>--</td><td>1.000000e+02</td><td>1.000000e+02</td><td>1.000000e+03</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_DT(D,L_Disk)</td><td>1.305539e+19</td><td>--</td><td>--</td><td>--</td><td>1.305539e+19</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>1.000000e-01</td><td>--</td><td>--</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>1.000000e-01</td><td>--</td><td>--</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_in(D,L_Disk)</td><td>1.958309e+17</td><td>--</td><td>--</td><td>--</td><td>1.958309e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_out(D,R_BLR_in)</td><td>2.154140e+17</td><td>--</td><td>--</td><td>--</td><td>2.154140e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk(M)</td><td>4.261082e+45</td><td>--</td><td>--</td><td>--</td><td>4.261082e+45</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>2.991580e+04</td><td>2.991580e+04</td><td>1.307071e+04</td><td>--</td><td>3.018434e+04</td><td>1.000000e+04</td><td>1.000000e+05</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R(D,theta_open)</td><td>9.134742e+16</td><td>--</td><td>--</td><td>--</td><td>7.607512e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_H(M)</td><td>8.694991e+17</td><td>8.694991e+17</td><td>6.630953e+16</td><td>--</td><td>8.695426e+17</td><td>8.000000e+17</td><td>5.000000e+19</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>9.141046e-02</td><td>9.141046e-02</td><td>2.996898e-02</td><td>--</td><td>1.500000e-01</td><td>1.000000e-02</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>1.000000e+00</td><td>--</td><td>--</td><td>--</td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta</td><td>2.000000e+00</td><td>--</td><td>--</td><td>--</td><td>2.000000e+00</td><td>0.000000e+00</td><td>9.000000e+01</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>BulkFactor</td><td>1.842049e+01</td><td>1.842049e+01</td><td>3.696853e+01</td><td>--</td><td>2.000000e+01</td><td>1.000000e+01</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>5.930000e-01</td><td>--</td><td>--</td><td>--</td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta_open(M)</td><td>5.997353e+00</td><td>5.997353e+00</td><td>5.466490e+00</td><td>--</td><td>5.000000e+00</td><td>4.000000e+00</td><td>6.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table5594993216-572158').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5594993216-572158').dataTable({
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

    p=model_minimizer_lsb.plot_corr_matrix()



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_34_0.png


.. code:: ipython3

    %matplotlib inline
    fit_model.set_nu_grid(1E6,1E30,200)
    fit_model.eval()
    p2=fit_model.plot_model(sed_data=sed_data)
    p2.setlim(y_min=1E-14,y_max=1E-9,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_35_0.png


.. code:: ipython3

    from jetset.minimizer import ModelMinimizer
    model_minimizer_minuit=ModelMinimizer('minuit')
    #fit_model.freeze('jet_leptonic','theta_open')
    best_fit_minuit=model_minimizer_minuit.fit(fit_model,sed_data,3E10,1E29,fitname='EC-best-fit-minuit',repeat=2)


.. parsed-literal::

    filtering data in fit range = [3.000000e+10,1.000000e+29]
    data length 24
    ================================================================================
    
    *** start fit process ***
    ----- 
    fit run: 0



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    ====> simplex
    ====> migrad
    - best chisq=3.21399e+01
    
    fit run: 1
    - old chisq=3.21399e+01



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    ====> simplex
    ====> migrad
    - best chisq=3.16610e+01
    
    -------------------------------------------------------------------------
    Fit report
    
    Model: EC-best-fit-minuit



.. raw:: html

    <i>Table length=22</i>
    <table id="table5551469904-337574" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>3.893234e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>8.993352e+04</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>9.909851e+01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.734853e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>1.274681e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.503112e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>DT</td><td>K</td><td>5.893807e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>1.305539e+19</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.958309e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.154140e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>4.261082e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>Disk</td><td>K</td><td>2.778004e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>3.603492e+17</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>3.509696e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>5.312839e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>2.000000e+00</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.862158e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table5551469904-337574').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5551469904-337574').dataTable({
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
    calls=1021
    mesg=



.. raw:: html

    <table>
        <tr>
            <th colspan="5" style="text-align:center" title="Minimizer"> Migrad </th>
        </tr>
        <tr>
            <td colspan="2" style="text-align:left" title="Minimum value of function"> FCN = 31.66 </td>
            <td colspan="3" style="text-align:center" title="Total number of function and (optional) gradient evaluations"> Nfcn = 1021 </td>
        </tr>
        <tr>
            <td colspan="2" style="text-align:left" title="Estimated distance to minimum and goal"> EDM = 0.00019 (Goal: 0.0002) </td>
            <td colspan="3" style="text-align:center" title="Total run time of algorithms"> time = 31.7 sec </td>
        </tr>
        <tr>
            <td colspan="2" style="text-align:center;background-color:#92CCA6;color:black"> Valid Minimum </td>
            <td colspan="3" style="text-align:center;background-color:#FFF79A;color:black"> SOME Parameters at limit </td>
        </tr>
        <tr>
            <td colspan="2" style="text-align:center;background-color:#92CCA6;color:black"> Below EDM threshold (goal x 10) </td>
            <td colspan="3" style="text-align:center;background-color:#92CCA6;color:black"> Below call limit </td>
        </tr>
        <tr>
            <td style="text-align:center;background-color:#92CCA6;color:black"> Covariance </td>
            <td style="text-align:center;background-color:#92CCA6;color:black"> Hesse ok </td>
            <td style="text-align:center;background-color:#92CCA6;color:black" title="Is covariance matrix accurate?"> Accurate </td>
            <td style="text-align:center;background-color:#92CCA6;color:black" title="Is covariance matrix positive definite?"> Pos. def. </td>
            <td style="text-align:center;background-color:#92CCA6;color:black" title="Was positive definiteness enforced by Minuit?"> Not forced </td>
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
            <td> 3.89323 </td>
            <td> 0.00023 </td>
            <td>  </td>
            <td>  </td>
            <td> 2 </td>
            <td> 100 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 1 </th>
            <td> par_1 </td>
            <td> 89.93352e3 </td>
            <td> 0.00028e3 </td>
            <td>  </td>
            <td>  </td>
            <td> 1E+04 </td>
            <td> 1E+05 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 2 </th>
            <td> par_2 </td>
            <td> 99.0985 </td>
            <td> 0.0014 </td>
            <td>  </td>
            <td>  </td>
            <td> 0 </td>
            <td>  </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 3 </th>
            <td> par_3 </td>
            <td> 173.4853 </td>
            <td> 0.0017 </td>
            <td>  </td>
            <td>  </td>
            <td> 100 </td>
            <td> 500 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 4 </th>
            <td> par_4 </td>
            <td> 1.274681 </td>
            <td> 0.000007 </td>
            <td>  </td>
            <td>  </td>
            <td> 1 </td>
            <td> 2.5 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 5 </th>
            <td> par_5 </td>
            <td> 3.50311218 </td>
            <td> 0.00000004 </td>
            <td>  </td>
            <td>  </td>
            <td> 3 </td>
            <td> 4 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 6 </th>
            <td> par_6 </td>
            <td> 589 </td>
            <td> 14 </td>
            <td>  </td>
            <td>  </td>
            <td> 100 </td>
            <td> 1E+03 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 7 </th>
            <td> par_7 </td>
            <td> 27.8e3 </td>
            <td> 2.4e3 </td>
            <td>  </td>
            <td>  </td>
            <td> 1E+04 </td>
            <td> 1E+05 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 8 </th>
            <td> par_8 </td>
            <td> 3.50970e18 </td>
            <td> 0.00014e18 </td>
            <td>  </td>
            <td>  </td>
            <td> 8E+17 </td>
            <td> 5E+19 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 9 </th>
            <td> par_9 </td>
            <td> 53.1284e-3 </td>
            <td> 0.0027e-3 </td>
            <td>  </td>
            <td>  </td>
            <td> 0.01 </td>
            <td> 1 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 10 </th>
            <td> par_10 </td>
            <td> 10.0 </td>
            <td> 0.4 </td>
            <td>  </td>
            <td>  </td>
            <td> 10 </td>
            <td> 30 </td>
            <td>  </td>
        </tr>
        <tr>
            <th> 11 </th>
            <td> par_11 </td>
            <td> 5.862158 </td>
            <td> 0.000008 </td>
            <td>  </td>
            <td>  </td>
            <td> 4 </td>
            <td> 6 </td>
            <td>  </td>
        </tr>
    </table>


.. parsed-literal::

    dof=12
    chisq=31.661039, chisq/red=2.638420 null hypothesis sig=0.001561
    
    best fit pars



.. raw:: html

    <i>Table length=22</i>
    <table id="table5594999360-794602" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>3.893234e+00</td><td>3.893234e+00</td><td>2.336607e-04</td><td>--</td><td>6.675472e+00</td><td>2.000000e+00</td><td>1.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>8.993352e+04</td><td>8.993352e+04</td><td>2.756248e-01</td><td>--</td><td>1.581028e+04</td><td>1.000000e+04</td><td>1.000000e+05</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>9.909851e+01</td><td>9.909851e+01</td><td>1.397515e-03</td><td>--</td><td>5.485583e+02</td><td>0.000000e+00</td><td>--</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma_break</td><td>1.734853e+02</td><td>1.734853e+02</td><td>1.658699e-03</td><td>--</td><td>1.894755e+02</td><td>1.000000e+02</td><td>5.000000e+02</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>1.274681e+00</td><td>1.274681e+00</td><td>6.788564e-06</td><td>--</td><td>1.588509e+00</td><td>1.000000e+00</td><td>2.500000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p_1</td><td>3.503112e+00</td><td>3.503112e+00</td><td>3.994333e-08</td><td>--</td><td>3.466529e+00</td><td>3.000000e+00</td><td>4.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_DT</td><td>5.893807e+02</td><td>5.893807e+02</td><td>1.422676e+01</td><td>--</td><td>4.502140e+02</td><td>1.000000e+02</td><td>1.000000e+03</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R_DT(D,L_Disk)</td><td>1.305539e+19</td><td>--</td><td>--</td><td>--</td><td>1.305539e+19</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_DT</td><td>1.000000e-01</td><td>--</td><td>--</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>tau_BLR</td><td>1.000000e-01</td><td>--</td><td>--</td><td>--</td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_in(D,L_Disk)</td><td>1.958309e+17</td><td>--</td><td>--</td><td>--</td><td>1.958309e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>*R_BLR_out(D,R_BLR_in)</td><td>2.154140e+17</td><td>--</td><td>--</td><td>--</td><td>2.154140e+17</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>L_Disk(M)</td><td>4.261082e+45</td><td>--</td><td>--</td><td>--</td><td>4.261082e+45</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>T_Disk</td><td>2.778004e+04</td><td>2.778004e+04</td><td>2.427427e+03</td><td>--</td><td>2.991580e+04</td><td>1.000000e+04</td><td>1.000000e+05</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>*R(D,theta_open)</td><td>3.603492e+17</td><td>--</td><td>--</td><td>--</td><td>9.134742e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>R_H(M)</td><td>3.509696e+18</td><td>3.509696e+18</td><td>1.379374e+14</td><td>--</td><td>8.694991e+17</td><td>8.000000e+17</td><td>5.000000e+19</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>5.312839e-02</td><td>5.312839e-02</td><td>2.726890e-06</td><td>--</td><td>9.141046e-02</td><td>1.000000e-02</td><td>1.000000e+00</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>1.000000e+00</td><td>--</td><td>--</td><td>--</td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta</td><td>2.000000e+00</td><td>--</td><td>--</td><td>--</td><td>2.000000e+00</td><td>0.000000e+00</td><td>9.000000e+01</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>BulkFactor</td><td>1.000000e+01</td><td>1.000000e+01</td><td>3.890289e-01</td><td>--</td><td>1.842049e+01</td><td>1.000000e+01</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>5.930000e-01</td><td>--</td><td>--</td><td>--</td><td>5.930000e-01</td><td>0.000000e+00</td><td>--</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>theta_open(M)</td><td>5.862158e+00</td><td>5.862158e+00</td><td>8.085871e-06</td><td>--</td><td>5.997353e+00</td><td>4.000000e+00</td><td>6.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table5594999360-794602').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5594999360-794602').dataTable({
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

    p=model_minimizer_minuit.plot_corr_matrix()



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_37_0.png


.. code:: ipython3

    %matplotlib inline
    fit_model.set_nu_grid(1E6,1E30,500)
    fit_model.eval()
    p2=fit_model.plot_model(sed_data=sed_data)
    p2.setlim(y_min=1E-14,y_max=1E-9,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_38_0.png


.. code:: ipython3

    jet.energetic_report()



.. raw:: html

    <i>Table length=39</i>
    <table id="table5595844896-295890" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
    <tr><td>BulkLorentzFactor</td><td>jet-bulk-factor</td><td></td><td>1.000000e+01</td></tr>
    <tr><td>U_e</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>3.869527e-03</td></tr>
    <tr><td>U_p_cold</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.489725e-01</td></tr>
    <tr><td>U_B</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.123087e-04</td></tr>
    <tr><td>U_Synch</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>9.582487e-06</td></tr>
    <tr><td>U_Synch_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>9.614802e-01</td></tr>
    <tr><td>U_Disk</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>2.827817e-06</td></tr>
    <tr><td>U_BLR</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>2.188308e-07</td></tr>
    <tr><td>U_DT</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.331679e-03</td></tr>
    <tr><td>U_CMB</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Star</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Disk_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>9.199862e-04</td></tr>
    <tr><td>U_BLR_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>5.115778e-05</td></tr>
    <tr><td>U_DT_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>6.675126e-06</td></tr>
    <tr><td>U_CMB_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Star_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_seed_tot</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.344308e-03</td></tr>
    <tr><td>L_Sync_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>4.687658e+41</td></tr>
    <tr><td>L_SSC_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>1.279464e+41</td></tr>
    <tr><td>L_EC_Disk_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_BLR_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>1.203063e+39</td></tr>
    <tr><td>L_EC_DT_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>6.871948e+42</td></tr>
    <tr><td>L_EC_CMB_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_Star_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_Sync</td><td>jet Lum.</td><td>erg / s</td><td>1.166040e+43</td></tr>
    <tr><td>jet_L_SSC</td><td>jet Lum.</td><td>erg / s</td><td>3.182626e+42</td></tr>
    <tr><td>jet_L_EC_Disk</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_BLR</td><td>jet Lum.</td><td>erg / s</td><td>2.992582e+40</td></tr>
    <tr><td>jet_L_EC_Star</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_DT</td><td>jet Lum.</td><td>erg / s</td><td>1.709376e+44</td></tr>
    <tr><td>jet_L_EC_CMB</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_pp_gamma</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_rad</td><td>jet Lum.</td><td>erg / s</td><td>1.858105e+44</td></tr>
    <tr><td>jet_L_kin</td><td>jet Lum.</td><td>erg / s</td><td>1.859850e+47</td></tr>
    <tr><td>jet_L_tot</td><td>jet Lum.</td><td>erg / s</td><td>1.863075e+47</td></tr>
    <tr><td>jet_L_e</td><td>jet Lum.</td><td>erg / s</td><td>4.708616e+45</td></tr>
    <tr><td>jet_L_B</td><td>jet Lum.</td><td>erg / s</td><td>1.366623e+44</td></tr>
    <tr><td>jet_L_p_cold</td><td>jet Lum.</td><td>erg / s</td><td>1.812764e+47</td></tr>
    <tr><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td></tr>
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
        console.log("$('#table5595844896-295890').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5595844896-295890').dataTable({
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

.. note::
    Please, read the introduction and the caveat :ref:`for the Bayesian model fitting <bayesian_model_fitting>` to understand the MCMC sampler workflow.


.. code:: ipython3

    from jetset.mcmc import McmcSampler
    from jetset.minimizer import ModelMinimizer
    model_minimizer_minuit = ModelMinimizer.load_model('EC_model_minimizer_minuit.pkl')



.. parsed-literal::

    ===> setting C threads to 12
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H


.. code:: ipython3

    mcmc=McmcSampler(model_minimizer_minuit)


.. code:: ipython3

    labels=['N','B','BulkFactor','p_1','gamma_break']
    model_name='jet_leptonic'
    use_labels_dict={model_name:labels}
    mcmc.set_labels(use_labels_dict=use_labels_dict)

.. code:: ipython3

    mcmc.set_bounds(bound=5.0,bound_rel=True)


.. parsed-literal::

    par: N  best fit value:  109.64687686370732  mcmc bounds: [0, 657.881261182244]
    par: B  best fit value:  0.05025286584056864  mcmc bounds: [0.01, 0.30151719504341185]
    par: BulkFactor  best fit value:  10.000000488240653  mcmc bounds: [10, 30]
    par: p_1  best fit value:  3.5156351521693128  mcmc bounds: [3, 4]
    par: gamma_break  best fit value:  139.69119783722982  mcmc bounds: [100, 500]


.. code:: ipython3

    mcmc.run_sampler(nwalkers=20, burnin=50,steps=500,progress='notebook')


.. parsed-literal::

    ===> setting C threads to 12
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H
    mcmc run starting
    



.. parsed-literal::

      0%|          | 0/500 [00:00<?, ?it/s]


.. parsed-literal::

    mcmc run done, with 1 threads took 393.33 seconds


.. code:: ipython3

    print(mcmc.acceptance_fraction)


.. parsed-literal::

    0.5046


.. code:: ipython3

    mcmc.model.set_nu_grid(1E6,1E30,200)
    
    p=mcmc.plot_model(sed_data=sed_data,fit_range=[3E10, 1E27],size=100)
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_49_0.png


.. code:: ipython3

    p=mcmc.plot_model(sed_data=sed_data,fit_range=[3E10, 1E27],size=100,quantiles=[0.05,0.95])
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_50_0.png


To have a better rendering on the scatter plot, we redefine the plot
labels

.. code:: ipython3

    mcmc.labels




.. parsed-literal::

    ['N', 'B', 'BulkFactor', 'p_1', 'gamma_break']



.. code:: ipython3

    mcmc.set_plot_label('N',r'$N$')
    mcmc.set_plot_label('B',r'$B$')
    mcmc.set_plot_label('BulkFactor',r'$\Gamma$')
    mcmc.set_plot_label('p_1',r'$p_1$')
    mcmc.set_plot_label('gamma_break',r'$\gamma_{\rm break}$')

the code below lets you tuning the output

1) mpl.rcParams[‘figure.dpi’] if you increase it you get a better
   definition
2) title_fmt=“.2E” this is the format for python, 2 significant digits,
   scientific notation
3) title_kwargs=dict(fontsize=12) you can change the fontsize

.. code:: ipython3

    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 80
    f=mcmc.corner_plot(quantiles=(0.16, 0.5, 0.84),title_kwargs=dict(fontsize=12),title_fmt=".2E",use_math_text=True)




.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_55_0.png


.. code:: ipython3

    f=mcmc.plot_chain(log_plot=False)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_56_0.png


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


.. parsed-literal::

    ===> setting C threads to 12
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H
    ===> setting C threads to 12
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5
    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H


.. code:: ipython3

    ms.model.set_nu_grid(1E6,1E30,200)
    
    p=ms.plot_model(sed_data=sed_data,fit_range=[3E10, 1E27],size=100)
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_60_0.png


.. code:: ipython3

    p=ms.plot_model(sed_data=sed_data,fit_range=[3E10, 1E27],size=100,quantiles=[0.05,0.95])
    p.setlim(y_min=1E-13,x_min=1E6,x_max=2E28)



.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_61_0.png


.. code:: ipython3

    f=ms.plot_par('p_1',log_plot=False)




.. image:: Jet_example_model_fit_EC_files/Jet_example_model_fit_EC_62_0.png

