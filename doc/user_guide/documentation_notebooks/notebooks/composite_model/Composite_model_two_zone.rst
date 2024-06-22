.. _composite_model_two_zones:

Composite Models: two-zone jets
===============================

.. code:: ipython3

    from jetset.jet_model import Jet
    from jetset.plot_sedfit import PlotSED
    from jetset.model_manager import FitModel

Here I present how to implement a two-zone model using the composite models are handled by  the :class:`.FitModel` class.
In the next release, I will provide a flexible convenience method (and, hopefully, the possibility to have mutual components IC emission).

We setup a two-zone model for an EC scenario, using a conical jet.

Setting the two components
--------------------------

.. code:: ipython3

    from jetset.jet_model import Jet

.. code:: ipython3

    compact_jet=Jet(name='compact_jet',beaming_expr='bulk_theta')


.. parsed-literal::

    ===> setting C threads to 12


.. code:: ipython3

    compact_jet.add_EC_component(['EC_DT','EC_BLR'],disk_type='BB')

.. code:: ipython3

    #we make the jet conical
    compact_jet.make_conical_jet()


.. parsed-literal::

    adding par: R_H to  R
    adding par: theta_open to  R
    ==> par R is depending on ['R_H', 'theta_open'] according to expr:   R =
    np.tan(np.radians(theta_open))*R_H
    setting R_H to 5.715026151380672e+16


.. code:: ipython3

    #we set the proper dependencies for the EC fields
    compact_jet.set_EC_dependencies()


.. parsed-literal::

    adding par: L_Disk to  R_BLR_in
    ==> par R_BLR_in is depending on ['L_Disk'] according to expr:   R_BLR_in =
    3E17*(L_Disk/1E46)**0.5
    adding par: R_BLR_in to  R_BLR_out
    ==> par R_BLR_out is depending on ['R_BLR_in'] according to expr:   R_BLR_out =
    R_BLR_in*1.1
    adding par: L_Disk to  R_DT
    ==> par R_DT is depending on ['L_Disk'] according to expr:   R_DT =
    2E19*(L_Disk/1E46)**0.5


.. code:: ipython3

    compact_jet.show_model()


.. parsed-literal::

    
    --------------------------------------------------------------------------------
    model description: 
    --------------------------------------------------------------------------------
    type: Jet
    name: compact_jet  
    geometry: spherical  
    
    electrons distribution:
     type: plc  
     gamma energy grid size:  201
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
     normalization:  True
     log-values:  False
     ratio of cold protons to relativistic electrons: 1.000000e+00
    
    accretion disk:
     disk Type: BB
     L disk: 1.000000e+45 (erg/s)
     T disk: 1.000000e+05 (K)
     nu peak disk: 8.171810e+15 (Hz)
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  100
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sum, hidden: False
       name:Sync, state: self-abs
       name:Sync, hidden: False
       name:SSC, state: on
       name:SSC, hidden: False
       name:EC_DT, state: on
       name:EC_DT, hidden: False
       name:DT, state: on
       name:DT, hidden: False
       name:Disk, state: on
       name:Disk, hidden: False
       name:EC_BLR, state: on
       name:EC_BLR, hidden: False
    external fields transformation method: blob
    
    SED info:
     nu grid size jetkernel: 1000
     nu size: 500
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    --------------------------------------------------------------------------------



.. raw:: html

    <i>Table length=21</i>
    <table id="table4733481504-844576" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>compact_jet</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>5.715026e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>1.000000e-01</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>z_cosm</td><td>redshift</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>5.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>1.000000e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_Disk</td><td>Disk</td><td>K</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.000000e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table4733481504-844576').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table4733481504-844576').dataTable({
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


Now we clone the ``compact`` jet to generate the ``extended`` region

.. code:: ipython3

    extended_jet=compact_jet.clone()


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


Setting the composite model
---------------------------

.. code:: ipython3

    extended_jet.name='extended_jet'

.. code:: ipython3

    from jetset.model_manager import FitModel

.. code:: ipython3

    composite_jet=FitModel(nu_size=500, name='composite_jet')


.. parsed-literal::

    /Users/orion/miniforge3/envs/jetset/lib/python3.10/site-packages/jetset/model_manager.py:158: UserWarning: no cosmology defined, using FlatLambdaCDM(name="Planck13", H0=67.77 km / (Mpc s), Om0=0.30712, Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.   0.   0.06] eV, Ob0=0.048252)
      warnings.warn(m)


.. code:: ipython3

    composite_jet.add_component(compact_jet)

.. code:: ipython3

    composite_jet.add_component(extended_jet)

We set the functional expression for the model composition

.. code:: ipython3

    composite_jet.composite_expr='extended_jet+compact_jet'

.. code:: ipython3

    composite_jet.show_pars()



.. raw:: html

    <i>Table length=42</i>
    <table id="table4746647520-167472" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>compact_jet</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>5.715026e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>1.000000e-01</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>z_cosm</td><td>redshift</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>5.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>1.000000e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_Disk</td><td>Disk</td><td>K</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.000000e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>T_DT</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>5.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>1.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>2.000000e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>1.000000e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>T_Disk</td><td>Disk</td><td>K</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>5.715026e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>theta</td><td>jet-viewing-angle</td><td>deg</td><td>1.000000e-01</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>z_cosm</td><td>redshift</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.000000e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table4746647520-167472').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table4746647520-167472').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



Linking pars in the composite model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important:: since we want the two-zone to share the same jet, we link **all** the jet-related parameters in the two model. For a pure SSC scenario, of course, you will not have  the dusty tours, and the accretion disk. 

.. code:: ipython3

    linked_pars=['z_cosm','theta_open','theta','T_Disk','L_Disk','T_DT']
    for par in linked_pars:
        composite_jet.link_par(par_name=par,from_model='extended_jet',to_model='compact_jet')



.. parsed-literal::

    adding par: z_cosm to  z_cosm
    adding par: theta_open to  theta_open
    adding par: theta to  theta
    adding par: T_Disk to  T_Disk
    adding par: L_Disk to  L_Disk
    adding par: T_DT to  T_DT


.. important:: Spectral components which are the same (i.e. duplicated), as the DT an the Disk, must be hidden in the jet component were are linked.

.. code:: ipython3

    #we hide the duplicated components from DT and Disk
    composite_jet.extended_jet.spectral_components.DT.hidden=True
    composite_jet.extended_jet.spectral_components.Disk.hidden=True

.. code:: ipython3

    composite_jet.show_pars()



.. raw:: html

    <i>Table length=42</i>
    <table id="table4748446640-334411" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>compact_jet</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>5.715026e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>theta(M)</td><td>jet-viewing-angle</td><td>deg</td><td>1.000000e-01</td><td>0.000000e+00</td><td>9.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>z_cosm(M)</td><td>redshift</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_DT(M)</td><td>DT</td><td>K</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>6.324555e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>9.486833e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>1.043552e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>compact_jet</td><td>L_Disk(M)</td><td>Disk</td><td>erg / s</td><td>1.000000e+45</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>T_Disk(M)</td><td>Disk</td><td>K</td><td>1.000000e+05</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>compact_jet</td><td>theta_open(M)</td><td>user_defined</td><td>deg</td><td>5.000000e+00</td><td>1.000000e+00</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+02</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>T_DT(L,compact_jet)</td><td>DT</td><td>K</td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>*R_DT(D,L_Disk)</td><td>DT</td><td>cm</td><td>6.324555e+18</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>tau_DT</td><td>DT</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>tau_BLR</td><td>BLR</td><td></td><td>1.000000e-01</td><td>0.000000e+00</td><td>1.000000e+00</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>*R_BLR_in(D,L_Disk)</td><td>BLR</td><td>cm</td><td>9.486833e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>*R_BLR_out(D,R_BLR_in)</td><td>BLR</td><td>cm</td><td>1.043552e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>L_Disk(L,compact_jet)</td><td>Disk</td><td>erg / s</td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>T_Disk(L,compact_jet)</td><td>Disk</td><td>K</td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>*R(D,theta_open)</td><td>region_size</td><td>cm</td><td>5.000000e+15</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>R_H(M)</td><td>region_position</td><td>cm</td><td>5.715026e+16</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-01</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>theta(L,compact_jet)</td><td>jet-viewing-angle</td><td>deg</td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>BulkFactor</td><td>jet-bulk-factor</td><td>lorentz-factor*</td><td>1.000000e+01</td><td>1.000000e+00</td><td>1.000000e+05</td><td>False</td><td>False</td></tr>
    <tr><td>extended_jet</td><td>z_cosm(L,compact_jet)</td><td>redshift</td><td></td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>extended_jet</td><td>theta_open(L,compact_jet)</td><td>user_defined</td><td>deg</td><td>--</td><td>--</td><td>--</td><td>False</td><td>True</td></tr>
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
        console.log("$('#table4748446640-334411').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table4748446640-334411').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [4, 5, 6], type: "optionalnum"}]
        });
    });
    </script>



setting pars in the extended and compat components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now place the extended region at a larger scale compared to the compact one, i.e. we increase ``R_H``, and we scale also the particle density and the magnetic field ``B``.

.. note::  The magnetic field could be scaled automatically, setting in both the components (``compact_jet`` and ``extended_jet``) the same functional dependency as shown in :ref:`dep_pars`, for the depending B field case. 

.. code:: ipython3

    composite_jet.compact_jet.parameters.R_H.val=1E18

.. code:: ipython3

    composite_jet.extended_jet.parameters.R_H.val=5E19

.. code:: ipython3

    composite_jet.extended_jet.parameters.gamma_cut.val=1E3

.. code:: ipython3

    composite_jet.extended_jet.parameters.B.val=0.001
    composite_jet.extended_jet.parameters.N.val=1

.. code:: ipython3

    composite_jet.eval()

.. code:: ipython3

    composite_jet.extended_jet.spectral_components.DT.hidden=True
    composite_jet.extended_jet.spectral_components.Disk.hidden=True

.. code:: ipython3

    composite_jet.eval()
    p=composite_jet.plot_model(skip_components=False)
    p.setlim(y_min=1E-16)



.. image:: Composite_model_two_zone_files/Composite_model_two_zone_38_0.png

