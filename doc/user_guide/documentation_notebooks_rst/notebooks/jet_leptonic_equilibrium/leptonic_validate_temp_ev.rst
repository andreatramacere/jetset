.. _leptonic_equilibrium_jet_validation_guide:

Validation of the leptonic equilibrium against the Fokker-Plank equation solution
=================================================================================

In this tutorial we validate the integral solution for the :math:`e^{\pm}` equilibrium, used for the leptonic jet against the Fokker-Plank equation solution implemented in the  :class:`.JetTimeEvol` class
(see the :ref:`leptonic_equilibrium_jet_guide`, for more details)

Jet leptonic
------------

.. code:: ipython3

    from jetset.jet_model import Jet
    from jetset.jet_emitters_factory import InjEmittersFactory
    
    import numpy as np
    import matplotlib.pyplot as plt


.. code:: ipython3

    import jetset
    print('tested with',jetset.__version__)


.. parsed-literal::

    tested with 1.4.0rc0


.. code:: ipython3

    
    q_inj = InjEmittersFactory().create_inj_emitters('pl')
    q_inj.parameters.gmin.val=100
    
    j = Jet(emitters_distribution=q_inj, emitters_type='electrons', verbose=False)
    j.parameters.T_esc_e_primaries.val=2
    j.parameters.R.val=1E16
    j.parameters.B.val=.2
    
    j.parameters.p.val=2.3
    j.parameters.L_inj.val=1E40
    j.parameters.gmin.val=1
    


.. code:: ipython3

    j.parameters


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=11</i>
    <table id="table6064732448-268950" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_leptonic</td><td>R</td><td>region_size</td><td>cm</td><td>1.000000e+16</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_leptonic</td><td>beam_obj</td><td>beaming</td><td></td><td>1.000000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>z_cosm</td><td>redshift</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.300000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>T_esc_e_primaries</td><td>escape_time</td><td>R / c</td><td>2.000000e+00</td><td>1.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>jet_leptonic</td><td>L_inj</td><td>L_inj</td><td>erg / s</td><td>1.000000e+40</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
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
        datatables: 'https://cdn.datatables.net/2.1.8/js/dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table6064732448-268950').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table6064732448-268950').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>





.. parsed-literal::

    None



.. code:: ipython3

    j.eval()
    j.plot_model()




.. parsed-literal::

    <jetset.plot_sedfit.PlotSED at 0x1696882f0>




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_8_1.png


.. code:: ipython3

    %matplotlib inline
    p=j.emitters_distribution.plot()




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_9_0.png


.. code:: ipython3

    j.save_model('leptonic_eq.pkl')

.. code:: ipython3

    j=Jet.load_model('leptonic_eq.pkl')

.. code:: ipython3

    j.eval()
    j.plot_model()




.. parsed-literal::

    <jetset.plot_sedfit.PlotSED at 0x16a8a88c0>




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_12_1.png


setting up the JetTimeEvol model
--------------------------------

.. code:: ipython3

    gamma_e_evovled=np.copy(j.emitters_distribution.gamma_e)
    n_gamma_e_evovled=np.copy(j.emitters_distribution.n_gamma_e)
    gamma_e_inj=np.copy(j.emitters_distribution.gamma_e_inj)
    n_gamma_e_inj=np.copy(j.emitters_distribution.n_gamma_e_inj)

.. code:: ipython3

    from jetset.jet_emitters_factory import EmittersFactory
    from jetset.jet_emitters import InjEmittersArrayDistribution
    q_inj=InjEmittersArrayDistribution(name='array_distr',emitters_type='electrons',gamma_array=gamma_e_inj,n_gamma_array=n_gamma_e_inj,normalize=False)

.. code:: ipython3

    q_inj.parameters


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=3</i>
    <table id="table6078571936-44060" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>Q</td><td>emitters_density</td><td>1 / (s cm3)</td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
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
        datatables: 'https://cdn.datatables.net/2.1.8/js/dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table6078571936-44060').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table6078571936-44060').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [3, 4, 5], type: "optionalnum"}]
        });
    });
    </script>





.. parsed-literal::

    None



.. code:: ipython3

    %matplotlib inline
    q_inj.update()
    p=q_inj.plot()
    p.ax.plot(gamma_e_inj, n_gamma_e_inj,'.',ms=1.5)
    





.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x16bd51f70>]




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_17_1.png


.. code:: ipython3

    from jetset.jet_timedep import JetTimeEvol
    temp_ev=JetTimeEvol(jet_rad=j,Q_inj=q_inj,only_radiation=True,inplace=True)


.. parsed-literal::

    /Users/orion/miniforge3/envs/py3.12/lib/python3.12/site-packages/jetset/model_manager.py:158: UserWarning: no cosmology defined, using FlatLambdaCDM(name="Planck13", H0=67.77 km / (Mpc s), Om0=0.30712, Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.   0.   0.06] eV, Ob0=0.048252)
      warnings.warn(m)


we use the acc region with escape time equal to radiative region

.. code:: ipython3

    duration=1E7
    duration_acc=0
    T_SIZE=np.int32(2E4)
    
    temp_ev.parameters.duration.val=duration
    
    temp_ev.parameters.TStart_Inj.val=0
    temp_ev.parameters.TStop_Inj.val=duration
    temp_ev.parameters.T_esc_rad.val= 2
    
    
    temp_ev.parameters.Esc_Index_rad.val=0
    temp_ev.parameters.t_size.val=T_SIZE
    temp_ev.parameters.num_samples.val=500
    temp_ev.IC_cooling='on'
    temp_ev.parameters.L_inj.val=0
    
    temp_ev.parameters.gmin_grid.val=1.1
    temp_ev.parameters.gmax_grid.val=5E9
    temp_ev.parameters.gamma_grid_size.val=400
    
    temp_ev.init_TempEv()
    temp_ev.region_expansion='off'
    temp_ev.show_model()


.. parsed-literal::

    --------------------------------------------------------------------------------
    JetTimeEvol model description
    --------------------------------------------------------------------------------
     
    physical setup: 
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=12</i>
    <table id="table6105341232-784690" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>val</th><th>units</th><th>val*</th><th>units*</th><th>log</th></tr></thead>
    <tr><td>delta t</td><td>time</td><td>5.000000e+02</td><td>s</td><td>0.0014989622899999999</td><td>R/c</td><td>False</td></tr>
    <tr><td>log. sampling</td><td>time</td><td>0.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R/c</td><td>time</td><td>3.335641e+05</td><td>s</td><td>1.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>IC cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Sync cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Adiab. cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Reg. expansion</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc rad</td><td>time</td><td>6.671282e+05</td><td>s</td><td>2.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>R_rad rad start</td><td>region_position</td><td>1.000000e+16</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_H rad start</td><td>region_position</td><td>1.000000e+17</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T min. synch. cooling</td><td></td><td>3.869000e+00</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>L inj (electrons)</td><td>injected lum.</td><td>9.998146e+39</td><td>erg/s</td><td>None</td><td></td><td>False</td></tr>
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
        datatables: 'https://cdn.datatables.net/2.1.8/js/dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table6105341232-784690').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table6105341232-784690').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [], type: "optionalnum"}]
        });
    });
    </script>



.. parsed-literal::

    
    model parameters: 
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=17</i>
    <table id="table6104607936-626070" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_time_ev</td><td>duration</td><td>time_grid</td><td>s</td><td>1.000000e+07</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmin_grid</td><td>gamma_grid</td><td></td><td>1.100000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmax_grid</td><td>gamma_grid</td><td></td><td>5.000000e+09</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gamma_grid_size</td><td>gamma_grid</td><td></td><td>4.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Inj</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Inj</td><td>time_grid</td><td>s</td><td>1.000000e+07</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_rad</td><td>escape_time</td><td>(R/c)*</td><td>2.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_rad</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_rad_start</td><td>region_size</td><td>cm</td><td>1.000000e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_H_rad_start</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>m_B</td><td>magnetic_field_index</td><td></td><td>1.000000e+00</td><td>1.000000e+00</td><td>2.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_jet_exp</td><td>exp_start_time</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>beta_exp_R</td><td>beta_expansion</td><td>v/c*</td><td>1.000000e+00</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_rad</td><td>magnetic_field</td><td>G</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_size</td><td>time_grid</td><td></td><td>2.000000e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>num_samples</td><td>time_ev_output</td><td></td><td>5.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>L_inj</td><td>inj_luminosity</td><td>erg / s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
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
        datatables: 'https://cdn.datatables.net/2.1.8/js/dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table6104607936-626070').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table6104607936-626070').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    p=temp_ev.plot_pre_run_plot(dpi=100)




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_21_0.png


.. code:: ipython3

    p=temp_ev.plot_time_profile()



.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_22_0.png




we perform the evoltion only for injection and cooling, withou
acceleration region

.. code:: ipython3

    temp_ev.run(only_injection=True,cache_SEDs_acc=False,do_injection=True,cache_SEDs_rad=False)



.. parsed-literal::

    temporal evolution running



.. parsed-literal::

      0%|          | 0/20000 [00:00<?, ?it/s]


.. parsed-literal::

    temporal evolution completed


.. code:: ipython3

    p=temp_ev.plot_tempev_emitters(region='rad',loglog=False,energy_unit='gamma',pow=0,plot_Q_inj=True)
    p.ax.plot(gamma_e_evovled,n_gamma_e_evovled,'-',label='analytical solution',lw=4,color='gray',alpha=0.75,zorder=0)
    p.ax.axvline(j._blob.emitters.gamma_cooling_eq,ls='--',lw='1')
    p.ax.legend()
    p.setlim(y_min=1E-15,y_max=5000,x_min=2)




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_26_0.png


.. code:: ipython3

    m=n_gamma_e_evovled>0
    x_analytical=np.log10(gamma_e_evovled[m])
    y_analytical=np.log10(n_gamma_e_evovled[m])
    
    m=temp_ev.rad_region.time_sampled_emitters.n_gamma[-1]>0
    x_num=np.log10(temp_ev.rad_region.time_sampled_emitters.gamma[m])
    y_num=np.log10(temp_ev.rad_region.time_sampled_emitters.n_gamma[-1][m])
    
    y_analytical_interp = np.interp(x_num, x_analytical,y_analytical, left=np.nan, right=np.nan)
    
    m=~np.isnan(y_analytical_interp)
    m=np.logical_and(m,x_num>1)
    m=np.logical_and(m,x_num<5.9)
    
    y_analytical_interp=10**y_analytical_interp[m]
    x_out=x_num[m]
    y_num=10**y_num[m]

.. code:: ipython3

    
    d=np.fabs(y_analytical_interp-y_num)/y_num
    assert(all(d<0.2))

.. code:: ipython3

    plt.plot(x_out,d)




.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x16c1ac470>]




.. image:: leptonic_validate_temp_ev_files/leptonic_validate_temp_ev_29_1.png


