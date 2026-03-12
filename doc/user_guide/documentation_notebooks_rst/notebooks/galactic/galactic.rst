.. _galactic_guide:

Galactic object models
======================

Tow models have been added to handle galactic objects:
1) :class:`.GalacticBeamed` to hadle galatcit objects that are not beamed
2) :class:`.GalacticUnbeamed` to hadle galatcit objects that are  beamed

The redshift will be set 0, and not used, and a new parameter, ``DL_cm``
will be used in place.

The only difference, is that the ``GalacticBeamed`` has beaming set as
default to 10 and set free, whilst ``GalacticUnbeamed`` has beaming set
as default to 1 and frozen

.. code:: ipython3

    from jetset.jet_model import GalacticUnbeamed,GalacticBeamed
    from astropy import constants as const
    from astropy import units as u
    import matplotlib.pyplot as plt
    import numpy as np

.. code:: ipython3

    import jetset
    print('tested with',jetset.__version__)


.. parsed-literal::

    tested with 1.4.0rc0


Galactic Unbeamed leptonic example
----------------------------------

This is an example for a PWN.

**Please notice, that for a realistic model you should build a composite
model with Radio and Wind region, and define a custom emitting region. A
proper thermal external field will be provided in the next version. The
``CMB`` is the same used for jets**

.. code:: ipython3

    pwn=GalacticUnbeamed(emitters_distribution='bkn',verbose=False,emitters_type='electrons',distance=2*u.kpc,name='pwn')
    pwn.add_EC_component(['EC_CMB'])
    pwn.parameters.R.val=(2*u.pc).to('cm').value
    pwn.parameters.B.val=1E-4
    pwn.parameters.gamma_break.val=1E8
    pwn.parameters.gmax.val=8E9
    pwn.parameters.gmin.val=5E5
    pwn.parameters.p.val=3.2
    pwn.parameters.p_1.val=3.8
    pwn.parameters.N.val=1E-8
    pwn.show_model()


.. parsed-literal::

    
    --------------------------------------------------------------------------------
    model description: 
    --------------------------------------------------------------------------------
    type: GalacticUnbeamed
    name: pwn  
    geometry: spherical  
    
    electrons distribution:
     type: bkn  
     gamma energy grid size:  201
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
     normalization:  True
     log-values:  False
     ratio of cold protons to relativistic electrons: 1.000000e+00
    
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
       name:EC_CMB, state: on
       name:EC_CMB, hidden: False
    external fields transformation method: disk
    
    SED info:
     nu grid size jetkernel: 1000
     nu size: 500
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=10</i>
    <table id="table13582765136-627533" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>pwn</td><td>R</td><td>region_size</td><td>cm</td><td>6.171355e+18</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>pwn</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>5.000000e+05</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>8.000000e+09</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e-08</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+08</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>p</td><td>LE_spectral_slope</td><td></td><td>3.200000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.800000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>DL_cm</td><td>distance</td><td>cm</td><td>6.171355e+21</td><td>0.000000e+00</td><td>3.085678e+24</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table13582765136-627533').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13582765136-627533').dataTable({
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


.. code:: ipython3

    pwn.eval()
    p=pwn.plot_model(frame='src')



.. image:: galactic_files/galactic_11_0.png


.. code:: ipython3

    pwn.save_model('pwn.pkl')

.. code:: ipython3

    GalacticUnbeamed.load_model("pwn.pkl")


.. parsed-literal::

    
    --------------------------------------------------------------------------------
    model description: 
    --------------------------------------------------------------------------------
    type: GalacticUnbeamed
    name: pwn  
    geometry: spherical  
    
    electrons distribution:
     type: bkn  
     gamma energy grid size:  201
     gmin grid : 5.000000e+05
     gmax grid : 8.000000e+09
     normalization:  True
     log-values:  False
     ratio of cold protons to relativistic electrons: 1.000000e+00
    
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
       name:EC_CMB, state: on
       name:EC_CMB, hidden: False
    external fields transformation method: disk
    
    SED info:
     nu grid size jetkernel: 1000
     nu size: 500
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=12</i>
    <table id="table13608946096-915514" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>pwn</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>5.000000e+05</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>8.000000e+09</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e-08</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>gamma_break</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+08</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>p</td><td>LE_spectral_slope</td><td></td><td>3.200000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>p_1</td><td>HE_spectral_slope</td><td></td><td>3.800000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>R</td><td>region_size</td><td>cm</td><td>6.171355e+18</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>pwn</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-04</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>NH_cold_to_rel_e</td><td>cold_p_to_rel_e_ratio</td><td></td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>pwn</td><td>beam_obj</td><td>beaming</td><td></td><td>1.000000e+00</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>pwn</td><td>DL_cm</td><td>distance</td><td>cm</td><td>6.171355e+21</td><td>0.000000e+00</td><td>3.085678e+24</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table13608946096-915514').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13608946096-915514').dataTable({
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




.. parsed-literal::

    None



.. code:: ipython3

    pwn.energetic_report()


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=19</i>
    <table id="table13586066672-511938" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
    <tr><td>BulkLorentzFactor</td><td>jet-bulk-factor</td><td></td><td>1.000000e+00</td></tr>
    <tr><td>U_e</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>7.480812e-09</td></tr>
    <tr><td>U_p_cold</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.503276e-11</td></tr>
    <tr><td>U_B</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>3.978874e-10</td></tr>
    <tr><td>U_Synch</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.512730e-11</td></tr>
    <tr><td>U_Disk</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_BLR</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_DT</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_CMB</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>4.185297e-13</td></tr>
    <tr><td>U_Star</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_seed_tot</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.512730e-11</td></tr>
    <tr><td>L_Sync</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>2.170466e+38</td></tr>
    <tr><td>L_SSC</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>2.312959e+35</td></tr>
    <tr><td>L_EC_Disk</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_BLR</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_DT</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_CMB</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>2.139695e+35</td></tr>
    <tr><td>L_EC_Star</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
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
        datatables: 'https://cdn.datatables.net/2.1.8/js/dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table13586066672-511938').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13586066672-511938').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [3], type: "optionalnum"}]
        });
    });
    </script>



Galactic Unbeamed hadronic pp example
-------------------------------------

This is an example for a galactic unbeamed hadronic (pp) object

.. code:: ipython3

    gal_hadronic=GalacticUnbeamed(emitters_distribution='plc',verbose=False,emitters_type='protons',distance=2*u.kpc,name='gal_hadronic_unbeamed')
    gal_hadronic.parameters.R.val=1E18
    gal_hadronic.parameters.N.val=1000
    gal_hadronic.parameters.B.val=1E-3
    
    gal_hadronic.show_model()


.. parsed-literal::

    
    --------------------------------------------------------------------------------
    model description: 
    --------------------------------------------------------------------------------
    type: GalacticUnbeamed
    name: gal_hadronic_unbeamed  
    geometry: spherical  
    
    protons distribution:
     type: plc  
     gamma energy grid size:  201
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
     normalization:  True
     log-values:  False
    
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
       name:PP_gamma, state: on
       name:PP_gamma, hidden: False
       name:PP_neutrino_tot, state: on
       name:PP_neutrino_tot, hidden: False
       name:PP_neutrino_mu, state: on
       name:PP_neutrino_mu, hidden: False
       name:PP_neutrino_e, state: on
       name:PP_neutrino_e, hidden: False
       name:Bremss_ep, state: on
       name:Bremss_ep, hidden: False
    external fields transformation method: disk
    
    SED info:
     nu grid size jetkernel: 1000
     nu size: 500
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=10</i>
    <table id="table13585841904-155232" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>gal_hadronic_unbeamed</td><td>R</td><td>region_size</td><td>cm</td><td>1.000000e+18</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>T_esc_e_secondaries</td><td>escape_time</td><td>R / c</td><td>1.000000e+00</td><td>1.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>NH_pp</td><td>target_density</td><td>1 / cm3</td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_unbeamed</td><td>DL_cm</td><td>distance</td><td>cm</td><td>6.171355e+21</td><td>0.000000e+00</td><td>3.085678e+24</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table13585841904-155232').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13585841904-155232').dataTable({
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


.. code:: ipython3

    gal_hadronic.eval()
    p=gal_hadronic.plot_model(frame='src')
    p.setlim(y_min=1E33)



.. image:: galactic_files/galactic_18_0.png


.. code:: ipython3

    gal_hadronic.energetic_report()


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=20</i>
    <table id="table13585380848-763929" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
    <tr><td>BulkLorentzFactor</td><td>jet-bulk-factor</td><td></td><td>1.000000e+00</td></tr>
    <tr><td>U_e</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>3.408013e-08</td></tr>
    <tr><td>U_B</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>3.978874e-08</td></tr>
    <tr><td>U_p</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>2.386056e+01</td></tr>
    <tr><td>U_p_target</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.503276e-03</td></tr>
    <tr><td>U_Synch</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.108748e-10</td></tr>
    <tr><td>U_Disk</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_BLR</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_DT</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_CMB</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Star</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_seed_tot</td><td>Energy dens.  rest. frame</td><td>erg / cm3</td><td>1.108748e-10</td></tr>
    <tr><td>L_Sync</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>4.176989e+37</td></tr>
    <tr><td>L_SSC</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>2.253964e+34</td></tr>
    <tr><td>L_EC_Disk</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_BLR</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_DT</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_CMB</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_Star</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_pp_gamma</td><td>Lum.  rest. frame.</td><td>erg / s</td><td>1.220405e+40</td></tr>
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
        console.log("$('#table13585380848-763929').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13585380848-763929').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [3], type: "optionalnum"}]
        });
    });
    </script>



Galactic Beamed hadronic pp example
-----------------------------------

.. code:: ipython3

    from jetset.jet_model import GalacticUnbeamed,GalacticBeamed
    from astropy import constants as const
    from astropy import units as u
    import matplotlib.pyplot as plt
    import numpy as np
    
    gal_hadronic=GalacticBeamed(emitters_distribution='plc',verbose=False,emitters_type='protons',distance=2*u.kpc,name='gal_hadronic_beamed')
    gal_hadronic.parameters.R.val=1E18
    gal_hadronic.parameters.N.val=1000
    gal_hadronic.parameters.B.val=1E-3
    
    gal_hadronic.show_model()


.. parsed-literal::

    
    --------------------------------------------------------------------------------
    model description: 
    --------------------------------------------------------------------------------
    type: GalacticBeamed
    name: gal_hadronic_beamed  
    geometry: spherical  
    
    protons distribution:
     type: plc  
     gamma energy grid size:  201
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
     normalization:  True
     log-values:  False
    
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
       name:PP_gamma, state: on
       name:PP_gamma, hidden: False
       name:PP_neutrino_tot, state: on
       name:PP_neutrino_tot, hidden: False
       name:PP_neutrino_mu, state: on
       name:PP_neutrino_mu, hidden: False
       name:PP_neutrino_e, state: on
       name:PP_neutrino_e, hidden: False
       name:Bremss_ep, state: on
       name:Bremss_ep, hidden: False
    external fields transformation method: blob
    
    SED info:
     nu grid size jetkernel: 1000
     nu size: 500
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
    --------------------------------------------------------------------------------


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=12</i>
    <table id="table13592185840-985482" class="table-striped table-bordered table-condensed">
    <thead><tr><th>model name</th><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
    <tr><td>gal_hadronic_beamed</td><td>R</td><td>region_size</td><td>cm</td><td>1.000000e+18</td><td>1.000000e+03</td><td>1.000000e+30</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>R_H</td><td>region_position</td><td>cm</td><td>1.000000e+17</td><td>0.000000e+00</td><td>--</td><td>False</td><td>True</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>B</td><td>magnetic_field</td><td>gauss</td><td>1.000000e-03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>T_esc_e_secondaries</td><td>escape_time</td><td>R / c</td><td>1.000000e+00</td><td>1.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>beam_obj</td><td>beaming</td><td></td><td>1.000000e+01</td><td>1.000000e-04</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.000000e+00</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1.000000e+06</td><td>1.000000e+00</td><td>1.000000e+15</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>N</td><td>emitters_density</td><td>1 / cm3</td><td>1.000000e+03</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>NH_pp</td><td>target_density</td><td>1 / cm3</td><td>1.000000e+00</td><td>0.000000e+00</td><td>--</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>gamma_cut</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>1.000000e+04</td><td>1.000000e+00</td><td>1.000000e+09</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>p</td><td>LE_spectral_slope</td><td></td><td>2.000000e+00</td><td>-1.000000e+01</td><td>1.000000e+01</td><td>False</td><td>False</td></tr>
    <tr><td>gal_hadronic_beamed</td><td>DL_cm</td><td>distance</td><td>cm</td><td>6.171355e+21</td><td>0.000000e+00</td><td>3.085678e+24</td><td>False</td><td>False</td></tr>
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
        console.log("$('#table13592185840-985482').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13592185840-985482').dataTable({
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


.. code:: ipython3

    gal_hadronic.eval()
    p=gal_hadronic.plot_model(frame='src')
    p.setlim(y_min=1E36)



.. image:: galactic_files/galactic_22_0.png



.. code:: ipython3

    gal_hadronic.energetic_report()


.. parsed-literal::

    WARNING: AstropyDeprecationWarning: 'classic' backend for show_in_notebook() is deprecated as of 6.1. Instead, use the supported backend 'ipydatagrid'. [astropy.table.table]



.. raw:: html

    <i>Table length=40</i>
    <table id="table13643051904-843031" class="table-striped table-bordered table-condensed">
    <thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
    <tr><td>BulkLorentzFactor</td><td>jet-bulk-factor</td><td></td><td>1.000000e+01</td></tr>
    <tr><td>U_e</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>3.408013e-08</td></tr>
    <tr><td>U_B</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>3.978874e-08</td></tr>
    <tr><td>U_p</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>2.386056e+01</td></tr>
    <tr><td>U_p_target</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.503276e-03</td></tr>
    <tr><td>U_Synch</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.108748e-10</td></tr>
    <tr><td>U_Synch_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>1.108748e-06</td></tr>
    <tr><td>U_Disk</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_BLR</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_DT</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_CMB</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Star</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Disk_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_BLR_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_DT_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_CMB_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_Star_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.000000e+00</td></tr>
    <tr><td>U_seed_tot</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>1.108748e-10</td></tr>
    <tr><td>L_Sync_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>4.176989e+37</td></tr>
    <tr><td>L_SSC_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>2.253964e+34</td></tr>
    <tr><td>L_EC_Disk_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_BLR_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_DT_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_CMB_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_EC_Star_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>L_pp_gamma_rf</td><td>Lum. blob rest. frame.</td><td>erg / s</td><td>1.220405e+40</td></tr>
    <tr><td>jet_L_Sync</td><td>jet Lum.</td><td>erg / s</td><td>1.039013e+39</td></tr>
    <tr><td>jet_L_SSC</td><td>jet Lum.</td><td>erg / s</td><td>5.606665e+35</td></tr>
    <tr><td>jet_L_EC_Disk</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_BLR</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_Star</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_DT</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_EC_CMB</td><td>jet Lum.</td><td>erg / s</td><td>0.000000e+00</td></tr>
    <tr><td>jet_L_pp_gamma</td><td>jet Lum.</td><td>erg / s</td><td>3.035720e+41</td></tr>
    <tr><td>jet_L_rad</td><td>jet Lum.</td><td>erg / s</td><td>3.046115e+41</td></tr>
    <tr><td>jet_L_kin</td><td>jet Lum.</td><td>erg / s</td><td>2.235985e+50</td></tr>
    <tr><td>jet_L_tot</td><td>jet Lum.</td><td>erg / s</td><td>2.235985e+50</td></tr>
    <tr><td>jet_L_e</td><td>jet Lum.</td><td>erg / s</td><td>3.193665e+41</td></tr>
    <tr><td>jet_L_B</td><td>jet Lum.</td><td>erg / s</td><td>3.728622e+41</td></tr>
    <tr><td>jet_L_p</td><td>jet Lum.</td><td>erg / s</td><td>2.235985e+50</td></tr>
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
        console.log("$('#table13643051904-843031').dataTable()");
    
    jQuery.extend( jQuery.fn.dataTableExt.oSort, {
        "optionalnum-asc": astropy_sort_num,
        "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
    });
    
        $('#table13643051904-843031').dataTable({
            order: [],
            pageLength: 100,
            lengthMenu: [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            pagingType: "full_numbers",
            columnDefs: [{targets: [3], type: "optionalnum"}]
        });
    });
    </script>



.. code:: ipython3

    gal_hadronic.save_model('gal_hadronic_beamed.pkl')

.. code:: ipython3

    gal_hadronic_new=GalacticBeamed.load_model('gal_hadronic_beamed.pkl')
