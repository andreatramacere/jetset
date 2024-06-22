.. _temp_ev_two_zone_cooling_acc:

Temporal evolution, two zones, cooling+acc
==========================================

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')

.. code:: ipython3

    import matplotlib.pyplot as plt
    import numpy as np

.. code:: ipython3

    import jetset
    print('tested on jetset',jetset.__version__)


.. parsed-literal::

    tested on jetset 1.3.0rc9


This is a very preliminary documentation for the temporal evolution
capabilities of jetset. Here we show how to create a decopuled
radiative+acceleration region, and how to evolve the system in order to
generate both particle spectra, SEDs, and lightcurves

To have full understanding of the analysis presented in this tutorial, it is advised to read the paper Tramacere et al. (2011) [Tramacere2011]_ for the understanding of stochastic and first order acceleration, and  Tramacere et al (2022) [Tramacere2022]_ for the coupling of radiative and acceleration regions.

definition of the injected particle distribution (``q_inj``), and of the
jet model for the radiative region

.. code:: ipython3

    from jetset.jet_emitters_factory import InjEmittersFactory
    from jetset.jet_model import Jet
    jet_model=Jet()
    q_inj=InjEmittersFactory().create_inj_emitters('pl',emitters_type='electrons',normalize=True)
    q_inj.parameters.gmin.val=9
    q_inj.parameters.gmax.val=10
    q_inj.parameters.p.val=0.5
    
    jet_model.parameters.beam_obj.val=30
    jet_model.parameters.B.val=0.2
    jet_model.parameters.z_cosm.val=0.03
    jet_model.parameters.R.val=5E15
    



.. parsed-literal::

    ===> setting C threads to 12


here we set some relevant parameters taht will be described in detail in
the next version of the documentation

.. code:: ipython3

    flare_duration=1.0E5
    duration=flare_duration*10
    t_D0=1.5E5
    t_A0=2.5E4
    T_esc_rad=1E60
    L_inj=5.0E39
    E_acc_max=4E60
    Delta_R_acc_ratio=0.1
    B_ratio=1.0
    T_SIZE=2E4
    NUM_SET=500
    Diff_Index=2.0
    Acc_Index=1.0

Here we instantiate the ``JetTimeEvol`` object, passing the radiative
region jet model, and the injected particle class.

.. code:: ipython3

    from jetset.jet_timedep import JetTimeEvol
    temp_ev_acc=JetTimeEvol(jet_rad=jet_model,Q_inj=q_inj,inplace=True)


.. parsed-literal::

    ===> setting C threads to 12
    ===> setting C threads to 12
    adding par: z_cosm to  z_cosm


.. Note:: 
**The IC cooling is switched off, as default, to make the process faster**.  to switch on the IC cooling `temp_ev_acc.IC_cooling='on'`

Now, we setup some relevant parameters

.. code:: ipython3

    temp_ev_acc.rad_region.jet.nu_min=1E8
    temp_ev_acc.acc_region.jet.nu_min=1E8
    T_SIZE=np.int32(T_SIZE)
    
    if Delta_R_acc_ratio is not None:
        temp_ev_acc.parameters.Delta_R_acc.val=temp_ev_acc.parameters.R_rad_start.val*Delta_R_acc_ratio
    
    T_esc_acc=t_A0/(temp_ev_acc.parameters.Delta_R_acc.val/3E10)*2
    
    
    
    temp_ev_acc.parameters.duration.val=duration
    temp_ev_acc.parameters.TStart_Acc.val=0
    temp_ev_acc.parameters.TStop_Acc.val=flare_duration
    temp_ev_acc.parameters.TStart_Inj.val=0
    temp_ev_acc.parameters.TStop_Inj.val=flare_duration
    temp_ev_acc.parameters.T_esc_acc.val=T_esc_acc
    temp_ev_acc.parameters.T_esc_rad.val=T_esc_rad
    temp_ev_acc.parameters.t_D0.val=t_D0
    temp_ev_acc.parameters.t_A0.val=t_A0
    temp_ev_acc.parameters.Esc_Index_acc.val=Diff_Index-2
    temp_ev_acc.parameters.Esc_Index_rad.val=0
    temp_ev_acc.parameters.Acc_Index.val=Acc_Index
    temp_ev_acc.parameters.Diff_Index.val=Diff_Index
    temp_ev_acc.parameters.t_size.val=T_SIZE
    temp_ev_acc.parameters.num_samples.val=NUM_SET
    temp_ev_acc.parameters.E_acc_max.val=E_acc_max
    temp_ev_acc.parameters.L_inj.val=L_inj
    
    
    temp_ev_acc.parameters.gmin_grid.val=1.0
    temp_ev_acc.parameters.gmax_grid.val=1E8
    temp_ev_acc.parameters.gamma_grid_size.val=1500
    
    temp_ev_acc.parameters.B_acc.val=temp_ev_acc.rad_region.jet.parameters.B.val*B_ratio
    temp_ev_acc.init_TempEv()
    temp_ev_acc.show_model()



.. parsed-literal::

    --------------------------------------------------------------------------------
    JetTimeEvol model description
    --------------------------------------------------------------------------------
     
    physical setup: 
    
    --------------------------------------------------------------------------------



.. raw:: html

    <i>Table length=29</i>
    <table id="table5785058896-369048" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>val</th><th>units</th><th>val*</th><th>units*</th><th>log</th></tr></thead>
    <tr><td>delta t</td><td>time</td><td>5.000000e+01</td><td>s</td><td>0.00029979245799999996</td><td>R/c</td><td>False</td></tr>
    <tr><td>log. sampling</td><td>time</td><td>0.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R/c</td><td>time</td><td>1.667820e+05</td><td>s</td><td>1.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>IC cooling</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Sync cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Adiab. cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Reg. expansion</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Diff coeff</td><td></td><td>6.666667e-06</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Acc coeff</td><td></td><td>4.000000e-05</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Diff index</td><td></td><td>2.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Acc index</td><td></td><td>1.000000e+00</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc acc</td><td>time</td><td>5.003461e+04</td><td>s</td><td>3.0</td><td>R_acc/c</td><td>False</td></tr>
    <tr><td>Eacc max</td><td>energy</td><td>4.000000e+60</td><td>erg</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc rad</td><td>time</td><td>1.667820e+65</td><td>s</td><td>1e+60</td><td>R/c</td><td>False</td></tr>
    <tr><td>Delta R acc</td><td>accelerator_width</td><td>5.000000e+14</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>B acc</td><td>magnetic field</td><td>2.000000e-01</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_rad rad start</td><td>region_position</td><td>5.000000e+15</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_H rad start</td><td>region_position</td><td>1.000000e+17</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T_A0=1/ACC_COEFF</td><td>time</td><td>2.500000e+04</td><td>s</td><td>0.149896229</td><td>R/c</td><td>False</td></tr>
    <tr><td>T_D0=1/DIFF_COEFF</td><td>time</td><td>1.500000e+05</td><td>s</td><td>0.899377374</td><td>R/c</td><td>False</td></tr>
    <tr><td>T_DA0=1/(2*DIFF_COEFF)</td><td>time</td><td>7.500000e+04</td><td>s</td><td>0.449688687</td><td>R/c</td><td>False</td></tr>
    <tr><td>gamma Lambda Turb.  max</td><td></td><td>1.173358e+11</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma Lambda Coher. max</td><td></td><td>1.173358e+10</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma eq Syst. Acc (synch. cool)</td><td></td><td>7.832383e+05</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma eq Diff. Acc (synch. cool)</td><td></td><td>1.309535e+05</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T cooling(gamma_eq=gamma_eq_Diff)</td><td></td><td>1.477242e+05</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T cooling(gamma_eq=gamma_eq_Sys)</td><td></td><td>2.469874e+04</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T min. synch. cooling</td><td></td><td>1.934500e+02</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>L inj (electrons)</td><td>injected lum.</td><td>5.000000e+39</td><td>erg/s</td><td>None</td><td></td><td>False</td></tr>
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
        console.log("$('#table5785058896-369048').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5785058896-369048').dataTable({
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



.. raw:: html

    <i>Table length=30</i>
    <table id="table5784500064-583824" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_time_ev</td><td>duration</td><td>time_grid</td><td>s</td><td>1.000000e+06</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmin_grid</td><td>gamma_grid</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmax_grid</td><td>gamma_grid</td><td></td><td>1.000000e+08</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gamma_grid_size</td><td>gamma_grid</td><td></td><td>1.500000e+03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Acc</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Acc</td><td>time_grid</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Inj</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Inj</td><td>time_grid</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_acc</td><td>escape_time</td><td>(R_acc/c)*</td><td>3.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_acc</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_D0</td><td>acceleration_time</td><td>s</td><td>1.500000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_A0</td><td>acceleration_time</td><td>s</td><td>2.500000e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Diff_Index</td><td>fp_coeff_index</td><td>s</td><td>2.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Acc_Index</td><td>fp_coeff_index</td><td></td><td>1.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Delta_R_acc</td><td>accelerator_width</td><td>cm</td><td>5.000000e+14</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_acc</td><td>magnetic_field</td><td>G</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>E_acc_max</td><td>acc_energy</td><td>erg</td><td>4.000000e+60</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Lambda_max_Turb</td><td>turbulence_scale</td><td>cm</td><td>1.000000e+15</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Lambda_choer_Turb_factor</td><td>turbulence_scale</td><td>cm</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_rad</td><td>escape_time</td><td>(R/c)*</td><td>1.000000e+60</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_rad</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_rad_start</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_H_rad_start</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>m_B</td><td>magnetic_field_index</td><td></td><td>1.000000e+00</td><td>1.000000e+00</td><td>2.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_jet_exp</td><td>exp_start_time</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>beta_exp_R</td><td>beta_expansion</td><td>v/c*</td><td>1.000000e+00</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_rad</td><td>magnetic_field</td><td>G</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_size</td><td>time_grid</td><td></td><td>2.000000e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>num_samples</td><td>time_ev_output</td><td></td><td>5.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>L_inj</td><td>inj_luminosity</td><td>erg / s</td><td>5.000000e+39</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
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
        console.log("$('#table5784500064-583824').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5784500064-583824').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    temp_ev_acc.plot_time_profile()




.. parsed-literal::

    <jetset.plot_sedfit.PlotTempEvDiagram at 0x158be6b60>




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_16_1.png


.. note:: 
    - **we do not want to evolve the particle in the `jet_rad`, so we set `only_injection=True`, and we set `do_injection=True` to injet the particle defined by `q_inj`**    
    - setting `cache_SEDs_rad=True` will generate and cache all the SED at any time of the `NUM_SET`. **This will increase the computational time during the run. Anyhow, will speed up the computation of SEDs and light curves. Moreover, these SEDs will be saved in the model, and read if you will reload the model in the future**.
    - setting `cache_SEDs_acc=True` will generate and cache also the SEDs in the acceleration region.

.. code:: ipython3

    only_injection=True
    do_injection=True
    plot_fit_model=True
    plot_fit_distr=True
    plot_emitters=True
    plot_lcs=True
    delta_t_out=1000
    eval_cross_time=False
    rest_frame='obs'
    temp_ev_acc.run(only_injection=only_injection,
                    do_injection=do_injection,
                    cache_SEDs_acc=True, 
                    cache_SEDs_rad=True)


.. parsed-literal::

    temporal evolution running



.. parsed-literal::

      0%|          | 0/20000 [00:00<?, ?it/s]


.. parsed-literal::

    temporal evolution completed
    caching SED for each saved distribution: start



.. parsed-literal::

      0%|          | 0/500 [00:00<?, ?it/s]


.. parsed-literal::

    caching SED for each saved distribution: done
    caching SED for each saved distribution: start



.. parsed-literal::

      0%|          | 0/500 [00:00<?, ?it/s]


.. parsed-literal::

    caching SED for each saved distribution: done


Particle spectrum in the radiative region

.. code:: ipython3

    p=temp_ev_acc.plot_tempev_emitters(region='rad',loglog=False,energy_unit='gamma',pow=0)
    p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_A, ls='--')
    p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_DA, ls='--')
    p.setlim(x_max=1E7,x_min=1,y_min=1E-18,y_max=100)



.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_20_0.png


Particle spectrum in the acceleration region

.. code:: ipython3

    p=temp_ev_acc.plot_tempev_emitters(region='acc',loglog=False,energy_unit='gamma',pow=0)
    p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_A, ls='--')
    p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_DA, ls='--')
    p.setlim(x_max=1E7,x_min=1,y_min=1E-30,y_max=100)




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_22_0.png


SEDs in the acceleration region

.. code:: ipython3

    p=temp_ev_acc.plot_tempev_model(region='rad',sed_data=None, use_cached = True)
    p.setlim(y_min=1E-18,x_min=1E7)



.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_24_0.png


SEDs in the acceleration region

.. code:: ipython3

    p=temp_ev_acc.plot_tempev_model(region='acc',sed_data=None, use_cached = True)
    p.setlim(y_min=1E-18,x_min=1E7)



.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_26_0.png


We generate a lightcurve in the range nu1=2.4E22 Hz, nu2=7.2E25 Hz,
without the effect of the light crossing time, in the observer frame

.. code:: ipython3

    lg=temp_ev_acc.rad_region.make_lc(nu1=2.4E22,nu2=7.2E25,name='gamma',eval_cross_time=False,delta_t_out=100,use_cached=True,frame='obs')


.. code:: ipython3

    lg




.. raw:: html

    <div><i>Table length=344</i>
    <table id="table5813176256" class="table-striped table-bordered table-condensed">
    <thead><tr><th>time</th><th>flux</th><th>R_blob</th><th>t_blob</th></tr></thead>
    <thead><tr><th>s</th><th>erg / (s cm2)</th><th>cm</th><th>s</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>0.0</td><td>0.0</td><td>5000000000000000.0</td><td>0.0</td></tr>
    <tr><td>100.0</td><td>0.0</td><td>5000000000000000.0</td><td>2912.6213592233007</td></tr>
    <tr><td>200.0</td><td>0.0</td><td>5000000000000000.0</td><td>5825.242718446601</td></tr>
    <tr><td>300.0</td><td>7.684633148016865e-86</td><td>5000000000000000.0</td><td>8737.864077669903</td></tr>
    <tr><td>400.0</td><td>3.183380180128441e-75</td><td>5000000000000000.0</td><td>11650.485436893203</td></tr>
    <tr><td>500.0</td><td>7.931984536411245e-61</td><td>5000000000000000.0</td><td>14563.106796116504</td></tr>
    <tr><td>600.0</td><td>6.94642582804193e-55</td><td>5000000000000000.0</td><td>17475.728155339806</td></tr>
    <tr><td>700.0</td><td>7.506058875474857e-47</td><td>5000000000000000.0</td><td>20388.349514563106</td></tr>
    <tr><td>800.0</td><td>6.3996495618931105e-43</td><td>5000000000000000.0</td><td>23300.970873786406</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>33400.0</td><td>1.4426695701183453e-10</td><td>5000000000000000.0</td><td>972815.5339805826</td></tr>
    <tr><td>33500.0</td><td>1.4339732703617928e-10</td><td>5000000000000000.0</td><td>975728.1553398059</td></tr>
    <tr><td>33600.0</td><td>1.4253399065766648e-10</td><td>5000000000000000.0</td><td>978640.7766990291</td></tr>
    <tr><td>33700.0</td><td>1.4167638069528754e-10</td><td>5000000000000000.0</td><td>981553.3980582524</td></tr>
    <tr><td>33800.0</td><td>1.4082471832460238e-10</td><td>5000000000000000.0</td><td>984466.0194174757</td></tr>
    <tr><td>33900.0</td><td>1.399789291822232e-10</td><td>5000000000000000.0</td><td>987378.6407766991</td></tr>
    <tr><td>34000.0</td><td>1.3913875298994561e-10</td><td>5000000000000000.0</td><td>990291.2621359223</td></tr>
    <tr><td>34100.0</td><td>1.3830459175970436e-10</td><td>5000000000000000.0</td><td>993203.8834951456</td></tr>
    <tr><td>34200.0</td><td>1.3747571962196792e-10</td><td>5000000000000000.0</td><td>996116.5048543689</td></tr>
    <tr><td>34300.0</td><td>1.3665298229570238e-10</td><td>5000000000000000.0</td><td>999029.1262135921</td></tr>
    </table></div>



.. code:: ipython3

    plt.plot(lg['time'],lg['flux'])
    plt.xlabel('time (%s)'%lg['time'].unit)
    plt.ylabel('flux (%s)'%lg['flux'].unit)




.. parsed-literal::

    Text(0, 0.5, 'flux (erg / (s cm2))')




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_30_1.png


We generate a lightcurve in the range nu1=2.4E22 Hz, nu2=7.2E25 Hz, with
the effect of the light crossing time, in the observer frame

.. code:: ipython3

    lg_cross=temp_ev_acc.rad_region.make_lc(nu1=2.4E22,nu2=7.2E25,name='gamma',eval_cross_time=True,delta_t_out=100,use_cached=True,frame='obs',cross_time_slices=100)


.. code:: ipython3

    plt.plot(lg['time'],lg['flux'])
    plt.plot(lg_cross['time'],lg_cross['flux'])
    
    plt.xlabel('time (%s)'%lg['time'].unit)
    plt.ylabel('flux (%s)'%lg['flux'].unit)




.. parsed-literal::

    Text(0, 0.5, 'flux (erg / (s cm2))')




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_33_1.png


We can save the model and reuse it later for plotting lightcurcves,
SEDs, and electron distributions

.. code:: ipython3

    temp_ev_acc.save_model('two_zone_rad_acc.pkl')

.. code:: ipython3

    temp_ev_acc_1=JetTimeEvol.load_model('two_zone_rad_acc.pkl')


.. parsed-literal::

    ===> setting C threads to 12
    ===> setting C threads to 12


.. code:: ipython3

    temp_ev_acc_1.show_model()


.. parsed-literal::

    --------------------------------------------------------------------------------
    JetTimeEvol model description
    --------------------------------------------------------------------------------
     
    physical setup: 
    
    --------------------------------------------------------------------------------



.. raw:: html

    <i>Table length=29</i>
    <table id="table4673812800-64039" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>val</th><th>units</th><th>val*</th><th>units*</th><th>log</th></tr></thead>
    <tr><td>delta t</td><td>time</td><td>5.000000e+01</td><td>s</td><td>0.00029979245799999996</td><td>R/c</td><td>False</td></tr>
    <tr><td>log. sampling</td><td>time</td><td>0.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R/c</td><td>time</td><td>1.667820e+05</td><td>s</td><td>1.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>IC cooling</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Sync cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Adiab. cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Reg. expansion</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Diff coeff</td><td></td><td>6.666667e-06</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Acc coeff</td><td></td><td>4.000000e-05</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Diff index</td><td></td><td>2.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Acc index</td><td></td><td>1.000000e+00</td><td>s-1</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc acc</td><td>time</td><td>5.003461e+04</td><td>s</td><td>3.0</td><td>R_acc/c</td><td>False</td></tr>
    <tr><td>Eacc max</td><td>energy</td><td>4.000000e+60</td><td>erg</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc rad</td><td>time</td><td>1.667820e+65</td><td>s</td><td>1e+60</td><td>R/c</td><td>False</td></tr>
    <tr><td>Delta R acc</td><td>accelerator_width</td><td>5.000000e+14</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>B acc</td><td>magnetic field</td><td>2.000000e-01</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_rad rad start</td><td>region_position</td><td>5.000000e+15</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_H rad start</td><td>region_position</td><td>1.000000e+17</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T_A0=1/ACC_COEFF</td><td>time</td><td>2.500000e+04</td><td>s</td><td>0.149896229</td><td>R/c</td><td>False</td></tr>
    <tr><td>T_D0=1/DIFF_COEFF</td><td>time</td><td>1.500000e+05</td><td>s</td><td>0.899377374</td><td>R/c</td><td>False</td></tr>
    <tr><td>T_DA0=1/(2*DIFF_COEFF)</td><td>time</td><td>7.500000e+04</td><td>s</td><td>0.449688687</td><td>R/c</td><td>False</td></tr>
    <tr><td>gamma Lambda Turb.  max</td><td></td><td>1.173358e+11</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma Lambda Coher. max</td><td></td><td>1.173358e+10</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma eq Syst. Acc (synch. cool)</td><td></td><td>7.832383e+05</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>gamma eq Diff. Acc (synch. cool)</td><td></td><td>1.309535e+05</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T cooling(gamma_eq=gamma_eq_Diff)</td><td></td><td>1.477242e+05</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T cooling(gamma_eq=gamma_eq_Sys)</td><td></td><td>2.469874e+04</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T min. synch. cooling</td><td></td><td>1.934500e+02</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>L inj (electrons)</td><td>injected lum.</td><td>5.000000e+39</td><td>erg/s</td><td>None</td><td></td><td>False</td></tr>
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
        console.log("$('#table4673812800-64039').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table4673812800-64039').dataTable({
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



.. raw:: html

    <i>Table length=30</i>
    <table id="table5794517568-970127" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_time_ev</td><td>duration</td><td>time_grid</td><td>s</td><td>1.000000e+06</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmin_grid</td><td>gamma_grid</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmax_grid</td><td>gamma_grid</td><td></td><td>1.000000e+08</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gamma_grid_size</td><td>gamma_grid</td><td></td><td>1.500000e+03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Acc</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Acc</td><td>time_grid</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Inj</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Inj</td><td>time_grid</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_acc</td><td>escape_time</td><td>(R_acc/c)*</td><td>3.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_acc</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_D0</td><td>acceleration_time</td><td>s</td><td>1.500000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_A0</td><td>acceleration_time</td><td>s</td><td>2.500000e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Diff_Index</td><td>fp_coeff_index</td><td>s</td><td>2.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Acc_Index</td><td>fp_coeff_index</td><td></td><td>1.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Delta_R_acc</td><td>accelerator_width</td><td>cm</td><td>5.000000e+14</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_acc</td><td>magnetic_field</td><td>G</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>E_acc_max</td><td>acc_energy</td><td>erg</td><td>4.000000e+60</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Lambda_max_Turb</td><td>turbulence_scale</td><td>cm</td><td>1.000000e+15</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Lambda_choer_Turb_factor</td><td>turbulence_scale</td><td>cm</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_rad</td><td>escape_time</td><td>(R/c)*</td><td>1.000000e+60</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_rad</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_rad_start</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_H_rad_start</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>m_B</td><td>magnetic_field_index</td><td></td><td>1.000000e+00</td><td>1.000000e+00</td><td>2.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_jet_exp</td><td>exp_start_time</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>beta_exp_R</td><td>beta_expansion</td><td>v/c*</td><td>1.000000e+00</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_rad</td><td>magnetic_field</td><td>G</td><td>2.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_size</td><td>time_grid</td><td></td><td>2.000000e+04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>num_samples</td><td>time_ev_output</td><td></td><td>5.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>L_inj</td><td>inj_luminosity</td><td>erg / s</td><td>5.000000e+39</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
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
        console.log("$('#table5794517568-970127').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5794517568-970127').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    p=temp_ev_acc_1.plot_tempev_model(region='rad',sed_data=None, use_cached = True)




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_38_0.png


.. code:: ipython3

    lx=temp_ev_acc_1.rad_region.make_lc(nu1=1E17,nu2=1E18,name='X',eval_cross_time=False,delta_t_out=100,use_cached=True,frame='obs')
    plt.plot(lx['time'],lx['flux'])
    plt.xlabel('time (%s)'%lg['time'].unit)
    plt.ylabel('flux (%s)'%lg['flux'].unit)




.. parsed-literal::

    Text(0, 0.5, 'flux (erg / (s cm2))')




.. image:: Temp_Ev_two_zones_acc_and_cooling_files/Temp_Ev_two_zones_acc_and_cooling_39_1.png


