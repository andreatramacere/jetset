.. _sherpa_plugin:

Example to use the sherpa plugin with the sherpa interface
==========================================================

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')
    
    import matplotlib.pylab as plt
    import jetset
    from jetset.test_data_helper import  test_SEDs
    from jetset.data_loader import ObsData,Data
    from jetset.plot_sedfit import PlotSED
    from jetset.test_data_helper import  test_SEDs

.. code:: ipython3

    test_SEDs




.. parsed-literal::

    ['/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_3C345.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk421_EBL_DEABS.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_ABS.ecsv',
     '/Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_DEABS.ecsv']



Loading data
------------

see the :ref:`data_format` user guide for further information about loading data 

.. code:: ipython3

    print(test_SEDs[2])
    data=Data.from_file(test_SEDs[2])



.. parsed-literal::

    /Users/orion/anaconda3/envs/jetset/lib/python3.8/site-packages/jetset/test_data/SEDs_data/SED_MW_Mrk501_EBL_ABS.ecsv


.. code:: ipython3

    %matplotlib inline
    sed_data=ObsData(data_table=data)
    sed_data.group_data(bin_width=0.2)
    
    sed_data.add_systematics(0.1,[10.**6,10.**29])
    p=sed_data.plot_sed()


.. parsed-literal::

    ================================================================================
    
    ***  binning data  ***
    ---> N bins= 90
    ---> bin_widht= 0.2
    ================================================================================
    



.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_7_1.png


.. code:: ipython3

    sed_data.save('Mrk_501.pkl')

phenomenological model constraining
-----------------------------------

see the :ref:`phenom_constr` user guide for further information about phenomenological constraining 

spectral indices
~~~~~~~~~~~~~~~~

.. code:: ipython3

    from jetset.sed_shaper import  SEDShape
    my_shape=SEDShape(sed_data)
    my_shape.eval_indices(minimizer='lsb',silent=True)
    p=my_shape.plot_indices()
    p.rescale(y_min=-15,y_max=-6)


.. parsed-literal::

    ================================================================================
    
    *** evaluating spectral indices for data ***
    ================================================================================
    



.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_12_1.png


sed shaper
~~~~~~~~~~

.. code:: ipython3

    mm,best_fit=my_shape.sync_fit(check_host_gal_template=True,
                      Ep_start=None,
                      minimizer='lsb',
                      silent=True,
                      fit_range=[10.,21.])


.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the synchrotron component ***
    ---> first blind fit run,  fit range: [10.0, 21.0]
    ---> class:  HSP
    
    ---> class:  HSP
    
    



.. raw:: html

    <i>Table length=6</i>
    <table id="table140682856225232-591923" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-6.411144e-02</td><td>-6.411144e-02</td><td>7.838965e-03</td><td>--</td><td>-4.778764e-02</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-1.751721e-03</td><td>-1.751721e-03</td><td>1.127030e-03</td><td>--</td><td>3.576201e-03</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>1.703747e+01</td><td>1.703747e+01</td><td>9.437354e-02</td><td>--</td><td>1.626870e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.030068e+01</td><td>-1.030068e+01</td><td>1.884114e-02</td><td>--</td><td>-1.025412e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>host_galaxy</td><td>nuFnu_p_host</td><td>-1.006557e+01</td><td>-1.006557e+01</td><td>5.462528e-02</td><td>--</td><td>-1.025412e+01</td><td>-1.225412e+01</td><td>-8.254123e+00</td><td>False</td></tr>
    <tr><td>host_galaxy</td><td>nu_scale</td><td>1.730764e-02</td><td>1.730764e-02</td><td>3.694887e-03</td><td>--</td><td>0.000000e+00</td><td>-5.000000e-01</td><td>5.000000e-01</td><td>False</td></tr>
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
        console.log("$('#table140682856225232-591923').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140682856225232-591923').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> sync       nu_p=+1.703747e+01 (err=+9.437354e-02)  nuFnu_p=-1.030068e+01 (err=+1.884114e-02) curv.=-6.411144e-02 (err=+7.838965e-03)
    ================================================================================
    


.. code:: ipython3

    my_shape.IC_fit(fit_range=[23.,29.],minimizer='minuit',silent=True)
    p=my_shape.plot_shape_fit()
    p.rescale(y_min=-15)


.. parsed-literal::

    ================================================================================
    
    *** Log-Polynomial fitting of the IC component ***
    ---> fit range: [23.0, 29.0]
    ---> LogCubic fit
    
    



.. raw:: html

    <i>Table length=4</i>
    <table id="table140682860988880-786788" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>val</th><th>bestfit val</th><th>err +</th><th>err -</th><th>start val</th><th>fit range min</th><th>fit range max</th><th>frozen</th></tr></thead>
    <tr><td>LogCubic</td><td>b</td><td>-1.565399e-01</td><td>-1.565399e-01</td><td>2.551779e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>c</td><td>-4.351917e-02</td><td>-4.351917e-02</td><td>2.032066e-02</td><td>--</td><td>-1.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Ep</td><td>2.529709e+01</td><td>2.529709e+01</td><td>1.817241e-01</td><td>--</td><td>2.536916e+01</td><td>0.000000e+00</td><td>3.000000e+01</td><td>False</td></tr>
    <tr><td>LogCubic</td><td>Sp</td><td>-1.058825e+01</td><td>-1.058825e+01</td><td>5.046950e-02</td><td>--</td><td>-1.000000e+01</td><td>-3.000000e+01</td><td>0.000000e+00</td><td>False</td></tr>
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
        console.log("$('#table140682860988880-786788').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140682860988880-786788').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [2, 3, 4, 5, 6, 7, 8], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    ---> IC         nu_p=+2.529709e+01 (err=+1.817241e-01)  nuFnu_p=-1.058825e+01 (err=+5.046950e-02) curv.=-1.565399e-01 (err=+2.551779e-02)
    ================================================================================
    



.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_15_3.png


Model constraining
~~~~~~~~~~~~~~~~~~

In this step we are not fitting the model, we are just obtaining the
phenomenological ``pre_fit`` model, that will be fitted in using minuit
ore least-square bound, as shown below

.. code:: ipython3

    from jetset.obs_constrain import ObsConstrain
    from jetset.model_manager import  FitModel
    sed_obspar=ObsConstrain(beaming=25,
                            B_range=[0.001,0.1],
                            distr_e='lppl',
                            t_var_sec=3*86400,
                            nu_cut_IR=1E12,
                            SEDShape=my_shape)
    
    
    prefit_jet=sed_obspar.constrain_SSC_model(electron_distribution_log_values=False,silent=True)
    prefit_jet.save_model('prefit_jet.pkl')


.. parsed-literal::

    ================================================================================
    
    ***  constrains parameters from observable ***
    



.. raw:: html

    <i>Table length=11</i>
    <table id="table140682861249536-286550" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>1.056955e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>5.050000e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>2.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>3.360000e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>4.703917e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>2.310708e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>7.087137e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma0_log_parab</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.045843e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>s</td><td>LE_spectral_slope</td><td></td><td>2.248787e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>r</td><td>spectral_curvature</td><td></td><td>3.205572e-01</td><td>-1.500000e+01</td><td>1.500000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table140682861249536-286550').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140682861249536-286550').dataTable({
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

    pl=prefit_jet.plot_model(sed_data=sed_data)
    pl.add_model_residual_plot(prefit_jet,sed_data)
    pl.rescale(y_min=-15,x_min=7,x_max=29)



.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_19_0.png


Model fitting with Sherpa
-------------------------

.. code:: ipython3

    from jetset.sherpa_plugin import JetsetSherpaModel


.. code:: ipython3

    from jetset.template_2Dmodel import EBLAbsorptionTemplate
    ebl_franceschini=EBLAbsorptionTemplate.from_name('Franceschini_2008')

.. code:: ipython3

    from jetset.jet_model import Jet
    prefit_jet=Jet.load_model('prefit_jet.pkl')




.. raw:: html

    <i>Table length=11</i>
    <table id="table140682858833760-535674" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>4.703917e+02</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>2.310708e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>7.087137e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gamma0_log_parab</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.045843e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>s</td><td>LE_spectral_slope</td><td></td><td>2.248787e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>r</td><td>spectral_curvature</td><td></td><td>3.205572e-01</td><td>-1.500000e+01</td><td>1.500000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>1.056955e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>5.050000e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td>lorentz-factor*</td><td>2.500000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>3.360000e-02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table140682858833760-535674').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table140682858833760-535674').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    sherpa_model_jet=JetsetSherpaModel(prefit_jet)
    sherpa_model_gal=JetsetSherpaModel(my_shape.host_gal)
    sherpa_model_ebl=JetsetSherpaModel(ebl_franceschini)
    



.. parsed-literal::

    jetset model name R renamed to  R_sh due to sherpa internal naming convention


.. code:: ipython3

    sherpa_model=(sherpa_model_jet+sherpa_model_gal)*sherpa_model_ebl

.. code:: ipython3

    sherpa_model




.. raw:: html

    <style>/*
    Copyright (C) 2020  Smithsonian Astrophysical Observatory
    
    
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License along
     with this program; if not, write to the Free Software Foundation, Inc.,
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    */
    
    :root {
      --sherpa-border-color: var(--jp-border-color2, #e0e0e0);
      --sherpa-background-color: var(--jp-layout-color0, white);
      --sherpa-background-color-row-even: var(--jp-layout-color1, white);
      --sherpa-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    
      /* https://medium.com/ge-design/iot-cool-gray-is-a-great-background-color-for-data-visualization-ebf18c318418 */
      --sherpa-background-color-dark1: #EBEFF2;
      --sherpa-background-color-dark2: #D8E0E5;
    }
    
    div.sherpa-text-fallback {
        display: none;
    }
    
    div.sherpa {
        display: block;
    }
    
    div.sherpa details summary {
        display: list-item;  /* needed for notebook, not lab */
        font-size: larger;
    }
    
    div.sherpa details div.datavals {
        display: grid;
        grid-template-columns: 1fr 3fr;
        column-gap: 0.5em;
    }
    
    div.sherpa div.dataname {
        font-weight: bold;
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa div.dataval { }
    
    div.sherpa div.datavals div:nth-child(4n + 1) ,
    div.sherpa div.datavals div:nth-child(4n + 2) {
        background: var(--sherpa-background-color-row-odd);
    }
    
    div.sherpa table.model tbody {
        border-bottom: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model tr.block {
        border-top: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd ,
    div.sherpa table.model th.model-even {
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd {
        background: var(--sherpa-background-color-dark1);
    }
    
    div.sherpa table.model th.model-even {
        background: var(--sherpa-background-color-dark2);
    }
    
    div.sherpa .failed {
        background: orange;
        font-size: large;
        padding: 1em;
    }
    </style><div class="sherpa-text-fallback">&lt;BinaryOpModel model instance &#x27;((jet_leptonic + host_galaxy) * Franceschini_2008)&#x27;&gt;</div><div hidden class="sherpa"><details open><summary>Model</summary><table class="model"><caption>Expression: (jet_leptonic + host_galaxy) * Franceschini_2008</caption><thead><tr><th>Component</th><th>Parameter</th><th>Thawed</th><th>Value</th><th>Min</th><th>Max</th><th>Units</th></tr></thead><tbody><tr><th class="model-odd" scope="rowgroup" rowspan=11>jet_leptonic</th><td>gmin</td><td><input disabled type="checkbox" checked></input></td><td>470.39174855643597</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>gmax</td><td><input disabled type="checkbox" checked></input></td><td>2310708.197406515</td><td>1.0</td><td>1000000000000000.0</td><td>lorentz-factor*</td></tr><tr><td>N</td><td><input disabled type="checkbox" checked></input></td><td>7.087136892517435</td><td>0.0</td><td>MAX</td><td>1 / cm3</td></tr><tr><td>gamma0_log_parab</td><td><input disabled type="checkbox" checked></input></td><td>10458.433237258416</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>s</td><td><input disabled type="checkbox" checked></input></td><td>2.2487867709713574</td><td>-10.0</td><td>10.0</td><td></td></tr><tr><td>r</td><td><input disabled type="checkbox" checked></input></td><td>0.32055721323507314</td><td>-15.0</td><td>15.0</td><td></td></tr><tr><td>R_sh</td><td><input disabled type="checkbox" checked></input></td><td>1.0569554329483122e+16</td><td>1000.0</td><td>1e+30</td><td>cm</td></tr><tr><td>R_H</td><td><input disabled type="checkbox" checked></input></td><td>1e+17</td><td>0.0</td><td>MAX</td><td>cm</td></tr><tr><td>B</td><td><input disabled type="checkbox" checked></input></td><td>0.0505</td><td>0.0</td><td>MAX</td><td>gauss</td></tr><tr><td>beam_obj</td><td><input disabled type="checkbox" checked></input></td><td>25.0</td><td>0.0001</td><td>MAX</td><td>lorentz-factor*</td></tr><tr><td>z_cosm</td><td><input disabled type="checkbox" checked></input></td><td>0.0336</td><td>0.0</td><td>MAX</td><td></td></tr><tr class="block"><th class="model-even" scope="rowgroup" rowspan=2>host_galaxy</th><td>nuFnu_p_host</td><td><input disabled type="checkbox" checked></input></td><td>-10.065573840795157</td><td>-12.25412262810351</td><td>-8.25412262810351</td><td>erg / (cm2 s)</td></tr><tr><td>nu_scale</td><td><input disabled type="checkbox" checked></input></td><td>0.01730764400081941</td><td>-0.5</td><td>0.5</td><td>Hz</td></tr><tr class="block"><th class="model-odd" scope="rowgroup" rowspan=1>Franceschini_2008</th><td>z_cosm</td><td><input disabled type="checkbox" checked></input></td><td>1.0</td><td>0.0</td><td>MAX</td><td></td></tr></tbody></table></details></div>



.. code:: ipython3

    sherpa_model_ebl.z_cosm  = sherpa_model_jet.z_cosm

.. code:: ipython3

    sherpa_model




.. raw:: html

    <style>/*
    Copyright (C) 2020  Smithsonian Astrophysical Observatory
    
    
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License along
     with this program; if not, write to the Free Software Foundation, Inc.,
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    */
    
    :root {
      --sherpa-border-color: var(--jp-border-color2, #e0e0e0);
      --sherpa-background-color: var(--jp-layout-color0, white);
      --sherpa-background-color-row-even: var(--jp-layout-color1, white);
      --sherpa-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    
      /* https://medium.com/ge-design/iot-cool-gray-is-a-great-background-color-for-data-visualization-ebf18c318418 */
      --sherpa-background-color-dark1: #EBEFF2;
      --sherpa-background-color-dark2: #D8E0E5;
    }
    
    div.sherpa-text-fallback {
        display: none;
    }
    
    div.sherpa {
        display: block;
    }
    
    div.sherpa details summary {
        display: list-item;  /* needed for notebook, not lab */
        font-size: larger;
    }
    
    div.sherpa details div.datavals {
        display: grid;
        grid-template-columns: 1fr 3fr;
        column-gap: 0.5em;
    }
    
    div.sherpa div.dataname {
        font-weight: bold;
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa div.dataval { }
    
    div.sherpa div.datavals div:nth-child(4n + 1) ,
    div.sherpa div.datavals div:nth-child(4n + 2) {
        background: var(--sherpa-background-color-row-odd);
    }
    
    div.sherpa table.model tbody {
        border-bottom: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model tr.block {
        border-top: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd ,
    div.sherpa table.model th.model-even {
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd {
        background: var(--sherpa-background-color-dark1);
    }
    
    div.sherpa table.model th.model-even {
        background: var(--sherpa-background-color-dark2);
    }
    
    div.sherpa .failed {
        background: orange;
        font-size: large;
        padding: 1em;
    }
    </style><div class="sherpa-text-fallback">&lt;BinaryOpModel model instance &#x27;((jet_leptonic + host_galaxy) * Franceschini_2008)&#x27;&gt;</div><div hidden class="sherpa"><details open><summary>Model</summary><table class="model"><caption>Expression: (jet_leptonic + host_galaxy) * Franceschini_2008</caption><thead><tr><th>Component</th><th>Parameter</th><th>Thawed</th><th>Value</th><th>Min</th><th>Max</th><th>Units</th></tr></thead><tbody><tr><th class="model-odd" scope="rowgroup" rowspan=11>jet_leptonic</th><td>gmin</td><td><input disabled type="checkbox" checked></input></td><td>470.39174855643597</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>gmax</td><td><input disabled type="checkbox" checked></input></td><td>2310708.197406515</td><td>1.0</td><td>1000000000000000.0</td><td>lorentz-factor*</td></tr><tr><td>N</td><td><input disabled type="checkbox" checked></input></td><td>7.087136892517435</td><td>0.0</td><td>MAX</td><td>1 / cm3</td></tr><tr><td>gamma0_log_parab</td><td><input disabled type="checkbox" checked></input></td><td>10458.433237258416</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>s</td><td><input disabled type="checkbox" checked></input></td><td>2.2487867709713574</td><td>-10.0</td><td>10.0</td><td></td></tr><tr><td>r</td><td><input disabled type="checkbox" checked></input></td><td>0.32055721323507314</td><td>-15.0</td><td>15.0</td><td></td></tr><tr><td>R_sh</td><td><input disabled type="checkbox" checked></input></td><td>1.0569554329483122e+16</td><td>1000.0</td><td>1e+30</td><td>cm</td></tr><tr><td>R_H</td><td><input disabled type="checkbox" checked></input></td><td>1e+17</td><td>0.0</td><td>MAX</td><td>cm</td></tr><tr><td>B</td><td><input disabled type="checkbox" checked></input></td><td>0.0505</td><td>0.0</td><td>MAX</td><td>gauss</td></tr><tr><td>beam_obj</td><td><input disabled type="checkbox" checked></input></td><td>25.0</td><td>0.0001</td><td>MAX</td><td>lorentz-factor*</td></tr><tr><td>z_cosm</td><td><input disabled type="checkbox" checked></input></td><td>0.0336</td><td>0.0</td><td>MAX</td><td></td></tr><tr class="block"><th class="model-even" scope="rowgroup" rowspan=2>host_galaxy</th><td>nuFnu_p_host</td><td><input disabled type="checkbox" checked></input></td><td>-10.065573840795157</td><td>-12.25412262810351</td><td>-8.25412262810351</td><td>erg / (cm2 s)</td></tr><tr><td>nu_scale</td><td><input disabled type="checkbox" checked></input></td><td>0.01730764400081941</td><td>-0.5</td><td>0.5</td><td>Hz</td></tr><tr class="block"><th class="model-odd" scope="rowgroup" rowspan=1>Franceschini_2008</th><td>z_cosm</td><td>linked</td><td>0.0336</td><td colspan=2>&#8656; jet_leptonic.z_cosm</td><td></td></tr></tbody></table></details></div>



.. code:: ipython3

    sherpa_model_jet.R_H.freeze()
    sherpa_model_jet.z_cosm.freeze()
    sherpa_model_gal.nu_scale.freeze()

.. code:: ipython3

    
    sherpa_model_jet.beam_obj.min = 5 
    sherpa_model_jet.beam_obj.max = 50.
    
    sherpa_model_jet.R_sh.min = 10**15. 
    sherpa_model_jet.R_sh.max = 10**17.5
    
    sherpa_model_jet.gmax.min = 1E5 
    sherpa_model_jet.gmax.max = 1E7
    


.. code:: ipython3

    sherpa_model




.. raw:: html

    <style>/*
    Copyright (C) 2020  Smithsonian Astrophysical Observatory
    
    
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License along
     with this program; if not, write to the Free Software Foundation, Inc.,
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    */
    
    :root {
      --sherpa-border-color: var(--jp-border-color2, #e0e0e0);
      --sherpa-background-color: var(--jp-layout-color0, white);
      --sherpa-background-color-row-even: var(--jp-layout-color1, white);
      --sherpa-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    
      /* https://medium.com/ge-design/iot-cool-gray-is-a-great-background-color-for-data-visualization-ebf18c318418 */
      --sherpa-background-color-dark1: #EBEFF2;
      --sherpa-background-color-dark2: #D8E0E5;
    }
    
    div.sherpa-text-fallback {
        display: none;
    }
    
    div.sherpa {
        display: block;
    }
    
    div.sherpa details summary {
        display: list-item;  /* needed for notebook, not lab */
        font-size: larger;
    }
    
    div.sherpa details div.datavals {
        display: grid;
        grid-template-columns: 1fr 3fr;
        column-gap: 0.5em;
    }
    
    div.sherpa div.dataname {
        font-weight: bold;
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa div.dataval { }
    
    div.sherpa div.datavals div:nth-child(4n + 1) ,
    div.sherpa div.datavals div:nth-child(4n + 2) {
        background: var(--sherpa-background-color-row-odd);
    }
    
    div.sherpa table.model tbody {
        border-bottom: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model tr.block {
        border-top: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd ,
    div.sherpa table.model th.model-even {
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd {
        background: var(--sherpa-background-color-dark1);
    }
    
    div.sherpa table.model th.model-even {
        background: var(--sherpa-background-color-dark2);
    }
    
    div.sherpa .failed {
        background: orange;
        font-size: large;
        padding: 1em;
    }
    </style><div class="sherpa-text-fallback">&lt;BinaryOpModel model instance &#x27;((jet_leptonic + host_galaxy) * Franceschini_2008)&#x27;&gt;</div><div hidden class="sherpa"><details open><summary>Model</summary><table class="model"><caption>Expression: (jet_leptonic + host_galaxy) * Franceschini_2008</caption><thead><tr><th>Component</th><th>Parameter</th><th>Thawed</th><th>Value</th><th>Min</th><th>Max</th><th>Units</th></tr></thead><tbody><tr><th class="model-odd" scope="rowgroup" rowspan=11>jet_leptonic</th><td>gmin</td><td><input disabled type="checkbox" checked></input></td><td>470.39174855643597</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>gmax</td><td><input disabled type="checkbox" checked></input></td><td>2310708.197406515</td><td>100000.0</td><td>10000000.0</td><td>lorentz-factor*</td></tr><tr><td>N</td><td><input disabled type="checkbox" checked></input></td><td>7.087136892517435</td><td>0.0</td><td>MAX</td><td>1 / cm3</td></tr><tr><td>gamma0_log_parab</td><td><input disabled type="checkbox" checked></input></td><td>10458.433237258416</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>s</td><td><input disabled type="checkbox" checked></input></td><td>2.2487867709713574</td><td>-10.0</td><td>10.0</td><td></td></tr><tr><td>r</td><td><input disabled type="checkbox" checked></input></td><td>0.32055721323507314</td><td>-15.0</td><td>15.0</td><td></td></tr><tr><td>R_sh</td><td><input disabled type="checkbox" checked></input></td><td>1.0569554329483122e+16</td><td>1000000000000000.0</td><td>3.1622776601683795e+17</td><td>cm</td></tr><tr><td>R_H</td><td><input disabled type="checkbox"></input></td><td>1e+17</td><td>0.0</td><td>MAX</td><td>cm</td></tr><tr><td>B</td><td><input disabled type="checkbox" checked></input></td><td>0.0505</td><td>0.0</td><td>MAX</td><td>gauss</td></tr><tr><td>beam_obj</td><td><input disabled type="checkbox" checked></input></td><td>25.0</td><td>5.0</td><td>50.0</td><td>lorentz-factor*</td></tr><tr><td>z_cosm</td><td><input disabled type="checkbox"></input></td><td>0.0336</td><td>0.0</td><td>MAX</td><td></td></tr><tr class="block"><th class="model-even" scope="rowgroup" rowspan=2>host_galaxy</th><td>nuFnu_p_host</td><td><input disabled type="checkbox" checked></input></td><td>-10.065573840795157</td><td>-12.25412262810351</td><td>-8.25412262810351</td><td>erg / (cm2 s)</td></tr><tr><td>nu_scale</td><td><input disabled type="checkbox"></input></td><td>0.01730764400081941</td><td>-0.5</td><td>0.5</td><td>Hz</td></tr><tr class="block"><th class="model-odd" scope="rowgroup" rowspan=1>Franceschini_2008</th><td>z_cosm</td><td>linked</td><td>0.0336</td><td colspan=2>&#8656; jet_leptonic.z_cosm</td><td></td></tr></tbody></table></details></div>



.. code:: ipython3

    from sherpa import data
    from sherpa.fit import Fit
    from sherpa.stats import Chi2
    from sherpa.optmethods import LevMar, NelderMead

.. code:: ipython3

    
    sherpa_data=data.Data1D("sed", sed_data.table['nu_data'], sed_data.table['nuFnu_data'], staterror=sed_data.table['dnuFnu_data'])


.. code:: ipython3

    fitter = Fit(sherpa_data, sherpa_model, stat=Chi2(), method=LevMar())
    fit_range=[1e11,1e29]
    
    sherpa_data.notice(fit_range[0], fit_range[1])


.. code:: ipython3

    results = fitter.fit()

.. code:: ipython3

    print("fit succeeded", results.succeeded)



.. parsed-literal::

    fit succeeded True


.. code:: ipython3

    results




.. raw:: html

    <style>/*
    Copyright (C) 2020  Smithsonian Astrophysical Observatory
    
    
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License along
     with this program; if not, write to the Free Software Foundation, Inc.,
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    */
    
    :root {
      --sherpa-border-color: var(--jp-border-color2, #e0e0e0);
      --sherpa-background-color: var(--jp-layout-color0, white);
      --sherpa-background-color-row-even: var(--jp-layout-color1, white);
      --sherpa-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    
      /* https://medium.com/ge-design/iot-cool-gray-is-a-great-background-color-for-data-visualization-ebf18c318418 */
      --sherpa-background-color-dark1: #EBEFF2;
      --sherpa-background-color-dark2: #D8E0E5;
    }
    
    div.sherpa-text-fallback {
        display: none;
    }
    
    div.sherpa {
        display: block;
    }
    
    div.sherpa details summary {
        display: list-item;  /* needed for notebook, not lab */
        font-size: larger;
    }
    
    div.sherpa details div.datavals {
        display: grid;
        grid-template-columns: 1fr 3fr;
        column-gap: 0.5em;
    }
    
    div.sherpa div.dataname {
        font-weight: bold;
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa div.dataval { }
    
    div.sherpa div.datavals div:nth-child(4n + 1) ,
    div.sherpa div.datavals div:nth-child(4n + 2) {
        background: var(--sherpa-background-color-row-odd);
    }
    
    div.sherpa table.model tbody {
        border-bottom: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model tr.block {
        border-top: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd ,
    div.sherpa table.model th.model-even {
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd {
        background: var(--sherpa-background-color-dark1);
    }
    
    div.sherpa table.model th.model-even {
        background: var(--sherpa-background-color-dark2);
    }
    
    div.sherpa .failed {
        background: orange;
        font-size: large;
        padding: 1em;
    }
    </style><div class="sherpa-text-fallback">&lt;Fit results instance&gt;</div><div hidden class="sherpa"><details open><summary>Fit parameters</summary><table class="fit"><thead><tr><th>Parameter</th><th>Best-fit value</th><th>Approximate error</th></tr></thead><tbody><tr><td>jet_leptonic.gmin</td><td>     221.318</td><td>&#177;      223.032</td></tr><tr><td>jet_leptonic.gmax</td><td> 2.12042e+06</td><td>&#177;            0</td></tr><tr><td>jet_leptonic.N</td><td>      9.7309</td><td>&#177;      9.52985</td></tr><tr><td>jet_leptonic.gamma0_log_parab</td><td>     6444.37</td><td>&#177;            0</td></tr><tr><td>jet_leptonic.s</td><td>     2.20226</td><td>&#177;    0.0932377</td></tr><tr><td>jet_leptonic.r</td><td>    0.225421</td><td>&#177;    0.0344871</td></tr><tr><td>jet_leptonic.R_sh</td><td> 1.47547e+16</td><td>&#177;            0</td></tr><tr><td>jet_leptonic.B</td><td>   0.0123178</td><td>&#177;   0.00255532</td></tr><tr><td>jet_leptonic.beam_obj</td><td>     41.9941</td><td>&#177;      5.52668</td></tr><tr><td>host_galaxy.nuFnu_p_host</td><td>    -10.0581</td><td>&#177;    0.0509935</td></tr></tbody></table></details><details><summary>Summary (10)</summary><div class="datavals"><div class="dataname">Method</div><div class="dataval">levmar</div><div class="dataname">Statistic</div><div class="dataval">chi2</div><div class="dataname">Final statistic</div><div class="dataval">8.25029</div><div class="dataname">Number of evaluations</div><div class="dataval">292</div><div class="dataname">Reduced statistic</div><div class="dataval">0.392871</div><div class="dataname">Probability (Q-value)</div><div class="dataval">0.993989</div><div class="dataname">Initial statistic</div><div class="dataval">124.614</div><div class="dataname">&#916; statistic</div><div class="dataval">116.364</div><div class="dataname">Number of data points</div><div class="dataval">31</div><div class="dataname">Degrees of freedom</div><div class="dataval">21</div></div></details></div>



.. code:: ipython3

    sherpa_model




.. raw:: html

    <style>/*
    Copyright (C) 2020  Smithsonian Astrophysical Observatory
    
    
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License along
     with this program; if not, write to the Free Software Foundation, Inc.,
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    
    */
    
    :root {
      --sherpa-border-color: var(--jp-border-color2, #e0e0e0);
      --sherpa-background-color: var(--jp-layout-color0, white);
      --sherpa-background-color-row-even: var(--jp-layout-color1, white);
      --sherpa-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    
      /* https://medium.com/ge-design/iot-cool-gray-is-a-great-background-color-for-data-visualization-ebf18c318418 */
      --sherpa-background-color-dark1: #EBEFF2;
      --sherpa-background-color-dark2: #D8E0E5;
    }
    
    div.sherpa-text-fallback {
        display: none;
    }
    
    div.sherpa {
        display: block;
    }
    
    div.sherpa details summary {
        display: list-item;  /* needed for notebook, not lab */
        font-size: larger;
    }
    
    div.sherpa details div.datavals {
        display: grid;
        grid-template-columns: 1fr 3fr;
        column-gap: 0.5em;
    }
    
    div.sherpa div.dataname {
        font-weight: bold;
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa div.dataval { }
    
    div.sherpa div.datavals div:nth-child(4n + 1) ,
    div.sherpa div.datavals div:nth-child(4n + 2) {
        background: var(--sherpa-background-color-row-odd);
    }
    
    div.sherpa table.model tbody {
        border-bottom: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model tr.block {
        border-top: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd ,
    div.sherpa table.model th.model-even {
        border-right: 1px solid var(--sherpa-border-color);
    }
    
    div.sherpa table.model th.model-odd {
        background: var(--sherpa-background-color-dark1);
    }
    
    div.sherpa table.model th.model-even {
        background: var(--sherpa-background-color-dark2);
    }
    
    div.sherpa .failed {
        background: orange;
        font-size: large;
        padding: 1em;
    }
    </style><div class="sherpa-text-fallback">&lt;BinaryOpModel model instance &#x27;((jet_leptonic + host_galaxy) * Franceschini_2008)&#x27;&gt;</div><div hidden class="sherpa"><details open><summary>Model</summary><table class="model"><caption>Expression: (jet_leptonic + host_galaxy) * Franceschini_2008</caption><thead><tr><th>Component</th><th>Parameter</th><th>Thawed</th><th>Value</th><th>Min</th><th>Max</th><th>Units</th></tr></thead><tbody><tr><th class="model-odd" scope="rowgroup" rowspan=11>jet_leptonic</th><td>gmin</td><td><input disabled type="checkbox" checked></input></td><td>221.31785661529102</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>gmax</td><td><input disabled type="checkbox" checked></input></td><td>2120423.887680081</td><td>100000.0</td><td>10000000.0</td><td>lorentz-factor*</td></tr><tr><td>N</td><td><input disabled type="checkbox" checked></input></td><td>9.7308970102023</td><td>0.0</td><td>MAX</td><td>1 / cm3</td></tr><tr><td>gamma0_log_parab</td><td><input disabled type="checkbox" checked></input></td><td>6444.369213670939</td><td>1.0</td><td>1000000000.0</td><td>lorentz-factor*</td></tr><tr><td>s</td><td><input disabled type="checkbox" checked></input></td><td>2.2022618914877032</td><td>-10.0</td><td>10.0</td><td></td></tr><tr><td>r</td><td><input disabled type="checkbox" checked></input></td><td>0.225421000374633</td><td>-15.0</td><td>15.0</td><td></td></tr><tr><td>R_sh</td><td><input disabled type="checkbox" checked></input></td><td>1.475471290427188e+16</td><td>1000000000000000.0</td><td>3.1622776601683795e+17</td><td>cm</td></tr><tr><td>R_H</td><td><input disabled type="checkbox"></input></td><td>1e+17</td><td>0.0</td><td>MAX</td><td>cm</td></tr><tr><td>B</td><td><input disabled type="checkbox" checked></input></td><td>0.012317793468494567</td><td>0.0</td><td>MAX</td><td>gauss</td></tr><tr><td>beam_obj</td><td><input disabled type="checkbox" checked></input></td><td>41.99405643997288</td><td>5.0</td><td>50.0</td><td>lorentz-factor*</td></tr><tr><td>z_cosm</td><td><input disabled type="checkbox"></input></td><td>0.0336</td><td>0.0</td><td>MAX</td><td></td></tr><tr class="block"><th class="model-even" scope="rowgroup" rowspan=2>host_galaxy</th><td>nuFnu_p_host</td><td><input disabled type="checkbox" checked></input></td><td>-10.058085094520706</td><td>-12.25412262810351</td><td>-8.25412262810351</td><td>erg / (cm2 s)</td></tr><tr><td>nu_scale</td><td><input disabled type="checkbox"></input></td><td>0.01730764400081941</td><td>-0.5</td><td>0.5</td><td>Hz</td></tr><tr class="block"><th class="model-odd" scope="rowgroup" rowspan=1>Franceschini_2008</th><td>z_cosm</td><td>linked</td><td>0.0336</td><td colspan=2>&#8656; jet_leptonic.z_cosm</td><td></td></tr></tbody></table></details></div>



.. code:: ipython3

    from jetset.sherpa_plugin import plot_sherpa_model


.. code:: ipython3

    p=plot_sherpa_model(sherpa_model_jet,label='SSC',line_style='--')
    p=plot_sherpa_model(sherpa_model_gal,plot_obj=p,label='host gal',line_style='--')
    p=plot_sherpa_model(sherpa_model=sherpa_model,plot_obj=p,sed_data=sed_data,fit_range=fit_range,add_res=True,label='(SSC+host gal)*ebl')
    




.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_40_0.png


You can access all the sherpa fetarues
https://sherpa.readthedocs.io/en/latest/fit/index.html

.. code:: ipython3

    from sherpa.plot import IntervalProjection
    iproj = IntervalProjection()
    iproj.prepare(fac=5, nloop=15)
    iproj.calc(fitter, sherpa_model_jet.s)
    iproj.plot()



.. image:: sherpa-plugin-sherpa-interface_files/sherpa-plugin-sherpa-interface_42_0.png


