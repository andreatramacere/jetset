.. _hadronic_pp_jet_validation_guide:

Validation of the pp equilibrium against the Fokker-Plank equation solution
===========================================================================

In this tutorial we validate the integral solution for the :math:`e^{\pm}` equilibrium, used for the pp jet against the Fokker-Plank equation solution implemented in the  :class:`.JetTimeEvol` class
(see the :ref:`hadronic_pp_jet_guide`, for more details on the pp hadronic jet model)

Jet pp
------

.. code:: ipython3

    from jetset.jet_model import Jet
    from jetset.jetkernel import jetkernel
    from astropy import constants as const
    from jetset.jet_emitters_factory import EmittersFactory
    from jetset.jet_emitters import InjEmittersArrayDistribution
    import numpy as np
    import matplotlib.pyplot as plt


.. code:: ipython3

    import jetset
    print('tested on jetset',jetset.__version__)


.. parsed-literal::

    tested on jetset 1.3.0rc7


.. code:: ipython3

    j=Jet(emitters_distribution='plc',verbose=False,emitters_type='protons')


.. parsed-literal::

    ===> setting C threads to 12


.. code:: ipython3

    
    j.parameters.z_cosm.val=z=0.001
    j.parameters.beam_obj.val=1
    j.parameters.gamma_cut.val=1000/(jetkernel.MPC2_TeV)
    j.parameters.NH_pp.val=1
    j.parameters.N.val=1
    j.parameters.p.val=2.0
    j.parameters.B.val=.5
    j.parameters.R.val=1E18
    j.parameters.gmin.val=1
    j.parameters.gmax.val=1E4
    j.set_emiss_lim(1E-60)
    j.set_IC_nu_size(100)
    j.gamma_grid_size=200
    j.eval()


.. code:: ipython3

    gmin=1.0/jetkernel.MPC2_TeV
    j.set_N_from_U_emitters(1.0, gmin=gmin)
    j.eval()
    
    #j.show_model()
    
    m=j.emitters_distribution.gamma_p>gmin
    print('U N(p) p>1 TeV=%e erg/cm-3'%(jetkernel.MPC2*np.trapz(j.emitters_distribution.n_gamma_p[m]*j.emitters_distribution.gamma_p[m],j.emitters_distribution.gamma_p[m])))


.. parsed-literal::

    U N(p) p>1 TeV=9.999992e-01 erg/cm-3


.. code:: ipython3

    j.energetic_report(verbose=False)

.. code:: ipython3

    %matplotlib inline
    j.emitters_distribution.plot()




.. parsed-literal::

    <jetset.plot_sedfit.PlotPdistr at 0x15357b7c0>




.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_10_1.png


.. code:: ipython3

    j.save_model('hadronic.pkl')

.. code:: ipython3

    from jetset.jet_model import Jet
    j=Jet.load_model('hadronic.pkl')


.. parsed-literal::

    ===> setting C threads to 12


setting up the JetTimeEvol model
--------------------------------

.. code:: ipython3

    gamma_sec_evovled=np.copy(j.emitters_distribution.gamma_e)
    n_gamma_sec_evovled=np.copy(j.emitters_distribution.n_gamma_e)
    gamma_sec_inj=np.copy(j.emitters_distribution.gamma_e_second_inj)
    n_gamma_sec_inj=np.copy(j.emitters_distribution.n_gamma_e_second_inj)

.. code:: ipython3

    from jetset.jet_emitters_factory import EmittersFactory
    from jetset.jet_emitters import InjEmittersArrayDistribution
    q_inj=InjEmittersArrayDistribution(name='array_distr',emitters_type='electrons',gamma_array=gamma_sec_inj,n_gamma_array=n_gamma_sec_inj,normalize=False)

.. code:: ipython3

    q_inj.parameters



.. raw:: html

    <i>Table length=3</i>
    <table id="table5713244896-349887" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.836150e+07</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
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
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table5713244896-349887').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5713244896-349887').dataTable({
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
    p=q_inj.plot()
    p.ax.plot(gamma_sec_inj, n_gamma_sec_inj,'.',ms=1.5)
    





.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x154a84b20>]




.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_17_1.png


.. code:: ipython3

    from jetset.jet_timedep import JetTimeEvol
    from jetset.jet_model import Jet
    
    temp_ev=JetTimeEvol(jet_rad=j,Q_inj=q_inj,only_radiation=True,inplace=True)


.. parsed-literal::

    ===> setting C threads to 12
    ===> setting C threads to 12


.. parsed-literal::

    /Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/model_manager.py:158: UserWarning: no cosmology defined, using FlatLambdaCDM(name="Planck13", H0=67.77 km / (Mpc s), Om0=0.30712, Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.   0.   0.06] eV, Ob0=0.048252)
      warnings.warn(m)


.. code:: ipython3

    temp_ev.Q_inj.parameters.Q.val




.. parsed-literal::

    1



we use the acc region with escape time equal to radiative region

.. code:: ipython3

    duration=5E9
    duration_acc=0
    T_SIZE=np.int32(2E6)
    
    temp_ev.parameters.duration.val=duration
    
    temp_ev.parameters.TStart_Inj.val=0
    temp_ev.parameters.TStop_Inj.val=duration
    temp_ev.parameters.T_esc_rad.val= 1
    
    
    temp_ev.parameters.Esc_Index_rad.val=0
    temp_ev.parameters.t_size.val=T_SIZE
    temp_ev.parameters.num_samples.val=500
    temp_ev.IC_cooling='off'
    temp_ev.parameters.L_inj.val=0
    
    temp_ev.parameters.gmin_grid.val=1.1
    temp_ev.parameters.gmax_grid.val=5E7
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



.. raw:: html

    <i>Table length=12</i>
    <table id="table5732819984-346763" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>par type</th><th>val</th><th>units</th><th>val*</th><th>units*</th><th>log</th></tr></thead>
    <tr><td>delta t</td><td>time</td><td>2.500000e+03</td><td>s</td><td>7.494811449999999e-05</td><td>R/c</td><td>False</td></tr>
    <tr><td>log. sampling</td><td>time</td><td>0.000000e+00</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R/c</td><td>time</td><td>3.335641e+07</td><td>s</td><td>1.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>IC cooling</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Sync cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Adiab. cooling</td><td></td><td>on</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Reg. expansion</td><td></td><td>off</td><td></td><td>None</td><td></td><td>False</td></tr>
    <tr><td>Tesc rad</td><td>time</td><td>3.335641e+07</td><td>s</td><td>1.0</td><td>R/c</td><td>False</td></tr>
    <tr><td>R_rad rad start</td><td>region_position</td><td>1.000000e+18</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>R_H rad start</td><td>region_position</td><td>1.000000e+17</td><td>cm</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>T min. synch. cooling</td><td></td><td>6.190400e+01</td><td>s</td><td>None</td><td></td><td>False</td></tr>
    <tr><td>L inj (electrons)</td><td>injected lum.</td><td>7.490567e+38</td><td>erg/s</td><td>None</td><td></td><td>False</td></tr>
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
        console.log("$('#table5732819984-346763').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5732819984-346763').dataTable({
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

    <i>Table length=17</i>
    <table id="table5732746336-555557" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>jet_time_ev</td><td>duration</td><td>time_grid</td><td>s</td><td>5.000000e+09</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmin_grid</td><td>gamma_grid</td><td></td><td>1.100000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gmax_grid</td><td>gamma_grid</td><td></td><td>5.000000e+07</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>gamma_grid_size</td><td>gamma_grid</td><td></td><td>4.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStart_Inj</td><td>time_grid</td><td>s</td><td>0.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>TStop_Inj</td><td>time_grid</td><td>s</td><td>5.000000e+09</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>T_esc_rad</td><td>escape_time</td><td>(R/c)*</td><td>1.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>Esc_Index_rad</td><td>fp_coeff_index</td><td></td><td>0.000000e+00</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_rad_start</td><td>region_size</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>R_H_rad_start</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>m_B</td><td>magnetic_field_index</td><td></td><td>1.000000e+00</td><td>1.000000e+00</td><td>2.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_jet_exp</td><td>exp_start_time</td><td>s</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>beta_exp_R</td><td>beta_expansion</td><td>v/c*</td><td>1.000000e+00</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>B_rad</td><td>magnetic_field</td><td>G</td><td>5.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>jet_time_ev</td><td>t_size</td><td>time_grid</td><td></td><td>2.000000e+06</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
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
        datatables: 'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table5732746336-555557').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table5732746336-555557').dataTable({
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




.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_22_0.png


.. code:: ipython3

    p=temp_ev.plot_time_profile()



.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_23_0.png


we perform the evoltion only for injection and cooling, withou
acceleration region

.. code:: ipython3

    temp_ev.run(only_injection=True,cache_SEDs_acc=False,do_injection=True,cache_SEDs_rad=False)



.. parsed-literal::

    temporal evolution running



.. parsed-literal::

      0%|          | 0/2000000 [00:00<?, ?it/s]


.. parsed-literal::

    temporal evolution completed


.. code:: ipython3

    p=temp_ev.plot_tempev_emitters(region='rad',loglog=False,energy_unit='gamma',pow=0,plot_Q_inj=True)
    p.ax.plot(gamma_sec_evovled,n_gamma_sec_evovled,'-',label='analytical solution',lw=4,color='gray',alpha=0.75,zorder=0)
    p.ax.legend()
    p.setlim(y_min=1E-30,y_max=1E-2)




.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_26_0.png


.. code:: ipython3

    m=n_gamma_sec_evovled>0
    x_analytical=np.log10(gamma_sec_evovled[m])
    y_analytical=np.log10(n_gamma_sec_evovled[m])
    
    m=temp_ev.rad_region.time_sampled_emitters.n_gamma[-1]>0
    x_num=np.log10(temp_ev.rad_region.time_sampled_emitters.gamma[m])
    y_num=np.log10(temp_ev.rad_region.time_sampled_emitters.n_gamma[-1][m])
    
    y_analytical_interp = np.interp(x_num, x_analytical,y_analytical, left=np.nan, right=np.nan)
    
    m=~np.isnan(y_analytical_interp)
    m=np.logical_and(m,x_num>0.25)
    m=np.logical_and(m,x_num<6)
    
    y_analytical_interp=10**y_analytical_interp[m]
    x_out=x_num[m]
    y_num=10**y_num[m]

.. code:: ipython3

    
    d=np.fabs(y_analytical_interp-y_num)/y_num
    assert(all(d<0.25))

.. code:: ipython3

    plt.plot(x_out,d)




.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x15580e920>]




.. image:: hadronic_validate_temp_ev_files/hadronic_validate_temp_ev_29_1.png

