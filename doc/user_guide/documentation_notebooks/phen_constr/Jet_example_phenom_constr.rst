.. _phenom_constr:

Phenomenological model constraining: application
================================================

.. figure:: ../slides/jetset_slides/jetset_slides.025.png
   :alt: image.png

   image.png

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')
    
    import matplotlib.pylab as plt
    from jetset.test_data_helper import  test_SEDs
    from jetset.data_loader import ObsData,Data
    from jetset.plot_sedfit import PlotSED
    from jetset.test_data_helper import  test_SEDs


.. code:: ipython3

    import jetset
    print('tested on jetset',jetset.__version__)


.. parsed-literal::

    tested on jetset 1.2.0rc6


.. code:: ipython3

    print(test_SEDs[1])
    data=Data.from_file(test_SEDs[2])



.. parsed-literal::

    /Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk421_EBL_DEABS.ecsv


.. code:: ipython3

    %matplotlib inline
    from jetset.cosmo_tools import Cosmo
    c=Cosmo()
    sed_data=ObsData(data_table=data,cosmo=c)
    sed_data.group_data(bin_width=0.2)
    sed_data.add_systematics(0.2,[10.**6,10.**29])



.. parsed-literal::

    ================================================================================
    
    ***  binning data  ***
    ---> N bins= 90
    ---> bin_widht= 0.2
    ================================================================================
    


.. code:: ipython3

    sed_data.save('Mrk_501.pkl')

.. code:: ipython3

    p=sed_data.plot_sed()



.. image:: Jet_example_phenom_constr_files/Jet_example_phenom_constr_8_0.png


.. code:: ipython3

    from jetset.sed_shaper import  SEDShape
    my_shape=SEDShape(sed_data)
    my_shape.eval_indices()
    p=my_shape.plot_indices()
    p.rescale(y_min=-15,y_max=-6)


.. parsed-literal::

    ================================================================================
    
    *** evaluating spectral indices for data ***
    ================================================================================
    



.. image:: Jet_example_phenom_constr_files/Jet_example_phenom_constr_9_1.png


.. code:: ipython3

    mm,best_fit=my_shape.sync_fit(check_host_gal_template=True,
                      Ep_start=None,
                      minimizer='minuit',
                      silent=True,
                      fit_range=[10,21])
    
    try:
        x,y,z,fig,ax=mm.minimizer.draw_contour('Ep','b')
    except:
        pass
    
    try:
        x,y,fig,ax=mm.minimizer.draw_profile('Ep')
    except:
        pass



.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the synchrotron component ***
    ---> first blind fit run,  fit range: [10, 21]
    ---> class:  HSP
    
    ---> class:  HSP
    
    



.. raw:: html

    <i>Table length=6</i>
    <table id="table140274936330752-668205" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-7.140660e-02</td><td>-7.140660e-02</td><td>1.337839e-02</td><td>--</td><td>-5.480219e-02</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-2.625704e-03</td><td>-2.625704e-03</td><td>2.018418e-03</td><td>--</td><td>3.829925e-03</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>1.694995e+01</td><td>1.694995e+01</td><td>1.504736e-01</td><td>--</td><td>1.603681e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.028879e+01</td><td>-1.028879e+01</td><td>3.653042e-02</td><td>--</td><td>-1.021025e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>host_galaxy</td><td>nuFnu_p_host</td><td>-1.007416e+01</td><td>-1.007416e+01</td><td>8.084766e-02</td><td>--</td><td>-1.021025e+01</td><td>-1.221025e+01</td><td>-8.210253e+00</td><td>False</td></tr>
    <tr><td>host_galaxy</td><td>nu_scale</td><td>-1.538032e-02</td><td>-1.538032e-02</td><td>3.055112e-05</td><td>--</td><td>0.000000e+00</td><td>-5.000000e-01</td><td>5.000000e-01</td><td>False</td></tr>
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
        console.log("$('#table140274936330752-668205').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140274936330752-668205').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> sync       nu_p=+1.694995e+01 (err=+1.504736e-01)  nuFnu_p=-1.028879e+01 (err=+3.653042e-02) curv.=-7.140660e-02 (err=+1.337839e-02)
    ================================================================================
    



.. image:: Jet_example_phenom_constr_files/Jet_example_phenom_constr_10_3.png


.. code:: ipython3

    help(mm.minimizer.minos_errors)


.. parsed-literal::

    Help on method minos_errors in module jetset.minimizer:
    
    minos_errors(par=None) method of jetset.minimizer.MinutiMinimizer instance
    


.. code:: ipython3

    my_shape.IC_fit(fit_range=[21,29],minimizer='lsb')
    p=my_shape.plot_shape_fit()
    p.rescale(y_min=-15,x_min=7,x_max=29)


.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the IC component ***
    ---> fit range: [21, 29]
    ---> LogCubic fit
    -------------------------------------------------------------------------
    Fit report
    
    Model: IC-shape-fit



.. raw:: html

    <i>Table length=4</i>
    <table id="table140274889755280-647348" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>curvature</td><td></td><td>-1.547080e-01</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>third-degree</td><td></td><td>-3.773118e-02</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>peak freq</td><td>Hz</td><td>2.526637e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>True</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>peak flux</td><td>erg / (cm2 s)</td><td>-1.057499e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>True</td><td>False</td></tr>
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
        console.log("$('#table140274889755280-647348').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140274889755280-647348').dataTable({
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
    calls=8
    mesg=



.. parsed-literal::

    'Both actual and predicted relative reductions in the sum of squares\n  are at most 0.000000 and the relative error between two consecutive iterates is at \n  most 0.000000'


.. parsed-literal::

    dof=9
    chisq=1.441196, chisq/red=0.160133 null hypothesis sig=0.997560
    
    stats without the UL
    dof  UL=9
    chisq=1.441196, chisq/red=0.160133 null hypothesis sig=0.997560
    
    
    best fit pars



.. raw:: html

    <i>Table length=4</i>
    <table id="table140274889803472-127894" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-1.547080e-01</td><td>-1.547080e-01</td><td>1.475869e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-3.773118e-02</td><td>-3.773118e-02</td><td>6.449603e-03</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>2.526637e+01</td><td>2.526637e+01</td><td>6.717304e-02</td><td>--</td><td>2.526894e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.057499e+01</td><td>-1.057499e+01</td><td>2.337071e-02</td><td>--</td><td>-1.000000e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table140274889803472-127894').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140274889803472-127894').dataTable({
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
    
    
    



.. raw:: html

    <i>Table length=4</i>
    <table id="table140274889804624-101986" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-1.547080e-01</td><td>-1.547080e-01</td><td>1.475869e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-3.773118e-02</td><td>-3.773118e-02</td><td>6.449603e-03</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>2.526637e+01</td><td>2.526637e+01</td><td>6.717304e-02</td><td>--</td><td>2.526894e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.057499e+01</td><td>-1.057499e+01</td><td>2.337071e-02</td><td>--</td><td>-1.000000e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table140274889804624-101986').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140274889804624-101986').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> IC         nu_p=+2.526637e+01 (err=+6.717304e-02)  nuFnu_p=-1.057499e+01 (err=+2.337071e-02) curv.=-1.547080e-01 (err=+1.475869e-02)
    ================================================================================
    



.. image:: Jet_example_phenom_constr_files/Jet_example_phenom_constr_12_9.png


.. code:: ipython3

    from jetset.obs_constrain import ObsConstrain
    from jetset.model_manager import  FitModel
    from jetset.minimizer import fit_SED
    sed_obspar=ObsConstrain(beaming=15,
                            B_range=[0.01,0.1],
                            distr_e='lppl',
                            t_var_sec=1*86400,
                            nu_cut_IR=5E10,
                            SEDShape=my_shape)
    
    
    jet=sed_obspar.constrain_SSC_model(electron_distribution_log_values=True,silent=False)


.. parsed-literal::

    ================================================================================
    
    ***  constrains parameters from observable ***
    
    ================================================================================
    
    ---> ***  emitting region parameters  ***
    
    ---> setting par type redshift, corresponding to par z_cosm
    
    ---> setting par type magnetic_field, corresponding to par B=5.500000e-02
    
    ---> setting par type region_size, corresponding to par R=3.759008e+16
    ---> completed True
    
    
    ---> *** electron distribution parameters ***
    ---> emitters distribution spectral type lp
    ---> emitters distribution name lppl
    
    ---> r elec. spec. curvature =3.570330e-01
    ---> setting par type curvature, corresponding to par r
    
    ---> s_radio_mm -0.47152657988709734 1.9430531597741947
    ---> s_X 3.269798782130266
    ---> s_Fermi 1.742749327549109
    ---> s_UV_X 2.745697034461969
    ---> s_Opt_UV -1.6299328530633286 4.259865706126657
    ---> s from synch log-log fit -1.0
    ---> s from (s_Fermi + s_UV)/2
    ---> power-law index s, class obj=HSP s chosen is 2.244223
    ---> setting par type LE_spectral_slope, corresponding to par s
    ---> task completed True
    
    ---> setting gamma_3p_Sync= 1.737092e+05, assuming B=5.500000e-02
    ---> task completed True
    
    ---> gamma_max=2.858471e+06 from nu_max_Sync= 2.413075e+19, using B=5.500000e-02
    ---> task completed True
    ---> setting par type high-energy-cut-off, corresponding to par gmax=6.456134e+00
    
    ---> setting par type low-energy-cut-off, corresponding to par gmin=2.114333e+00
    ---> task completed True
    
    ---> setting par type turn-over energy, corresponding to par gamma0_log_parab=4.181410e+00
    ---> task completed True
    ---> using gamma_3p_Sync= 173709.17000887552
    
    ---> nu_p_seed_blob=6.140587e+15
    ---> COMPTON FACTOR=8.632933e+00 19193.285984868435
    ---> determine gamma_3p_SSCc= 2.531574e+05
    ---> task completed True
    
    ---> setting par type turn-over energy, corresponding to par gamma0_log_parab=4.344977e+00
    ---> task completed True
    ---> using gamma_3p_SSC=2.531574e+05
    
    
    ---> setting par type emitters_density, corresponding to par N
    ---> to N=3.275365e+00
    ---> task completed (None, True)
    
    ---> setting B from nu_p_S to B=1.000000e+00
    ---> to B=1.000000e+00
    ---> setting B from best matching of nu_p_IC
    
         Best B=3.372112e-02
    ---> setting par type magnetic_field, corresponding to par B
    ---> task completed  True
    ---> best B found: 3.372112e-02
    
    ---> update pars for new B 
    ---> setting par type low-energy-cut-off, corresponding to par gmin
    ---> task completed True
    ---> set to 2.220564e+00
    
    ---> setting par type low-energy-cut-off, corresponding to par gamma0_log_parab
    ---> task completed True
    ---> task completed  True
    ---> using gamma_3p_Sync= 221846.7734643539
    ---> to 4.287640e+00
    
    ---> gamma_max=3.650600e+06 from nu_max_Sync= 2.413075e+19, using B=3.372112e-02
    ---> task completed True
    ---> setting par type high-energy-cut-off, corresponding to par gmax
    ---> set to 6.562364e+00
    
    ---> setting par type emitters_density, corresponding to par N
    ---> to N=7.006779e+00
    ---> task completed (None, True)
    
    ---> setting R from Compton Dominance (CD)
         Best R=3.942437e+16
    ---> setting par type region_size, corresponding to par R
    ---> set to 3.942437e+16
    ---> task completed True
    ---> updating setting par type emitters_density, corresponding to par N
    ---> set to 6.073564e+00
    ---> task completed (None, True)
    ---> t_var (days) 1.0487974107470437
    
    show pars



.. raw:: html

    <i>Table length=11</i>
    <table id="table140274940540912-436275" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>3.942437e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>3.372112e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>1.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>3.360000e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.220564e+00</td><td>0.000000e+00</td><td>9.000000e+00</td><td>True</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>6.562364e+00</td><td>0.000000e+00</td><td>1.500000e+01</td><td>True</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>6.073564e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma0_log_parab</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>4.287640e+00</td><td>0.000000e+00</td><td>9.000000e+00</td><td>True</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>s</td><td>LE_spectral_slope</td><td></td><td>2.244223e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>r</td><td>spectral_curvature</td><td></td><td>3.570330e-01</td><td>-1.500000e+01</td><td>1.500000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table140274940540912-436275').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140274940540912-436275').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    eval_model
    
    ================================================================================
    


.. code:: ipython3

    pl=jet.plot_model(sed_data=sed_data)
    pl.rescale(y_min=-15,x_min=7,x_max=29)
    jet.save_model('constrained_jet.pkl')



.. image:: Jet_example_phenom_constr_files/Jet_example_phenom_constr_14_0.png


