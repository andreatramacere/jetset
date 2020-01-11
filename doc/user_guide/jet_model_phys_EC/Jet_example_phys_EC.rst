.. _jet_physical_guide_EC:



External Compton
----------------


The external Compton implementation  gives you the possibility to use a double approach
 
* transformation of the external  fields to the blob rest frame :cite:`Dermer2000`

* transformation of the electron emitting distribution from the blob restframe to
  disk/BH restframe :cite:`Dermer95` :cite:`GKM01`

Regading the single external radiavite fiels 
 
* Implementation of Broad Line Region radiative field using the approach of :cite:`Donea2003` 

* Dusty torus implemented as a uniform BB field within `R_DT`

* accretion disk can be set a single BB or a multitemperature BB

* CMB 

Please read :ref:`jet_physical_guide_EC` if you skipped it.

.. figure:: jetset_EC_scheme.png
   :alt: EC scheme

   EC scheme

Broad Line Region
~~~~~~~~~~~~~~~~~

.. image::jetset_EC_scheme.png
  :width: 400
  :alt: EC scheme


.. code:: ipython3

    from jetset.jet_model import Jet
    my_jet=Jet(name='BLR example',electron_distribution='bkn',beaming_expr='bulk_theta')
    my_jet.add_EC_component(['EC_BLR','EC_Disk'])
    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: BLR example  
    
    electron distribution:
     type: bkn  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 1.000000e+06
     normalization  True
     log-values  False
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  50
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
       name:EC_BLR, state: on
       name:Disk, state: on
       name:EC_Disk, state: on
    external fields transformation method: blob
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
        name          par type           units             val         phys. bound. min  phys. bound. max   log  frozen
    ----------- ------------------- --------------- ------------------ ---------------- ------------------ ----- ------
              N    electron_density         1 / cm3              100.0                0               None False  False
           gmin  low-energy-cut-off lorentz-factor*                2.0                1       1000000000.0 False  False
           gmax high-energy-cut-off lorentz-factor*          1000000.0                1 1000000000000000.0 False  False
              p   LE_spectral_slope                                2.0            -10.0               10.0 False  False
            p_1   HE_spectral_slope                                3.0            -10.0               10.0 False  False
    gamma_break    turn-over-energy lorentz-factor*            10000.0                1       1000000000.0 False  False
              R         region_size              cm 5000000000000000.0           1000.0              1e+30 False  False
            R_H     region_position              cm              1e+17                0               None False   True
              B      magnetic_field               G                0.1                0               None False  False
          theta   jet-viewing-angle             deg                0.1                0               None False  False
     BulkFactor     jet-bulk-factor Lorentz-factor*               10.0              1.0               None False  False
         z_cosm            redshift                                0.1                0               None False  False
        tau_BLR                 BLR                                0.1                0                1.0 False  False
       R_BLR_in                 BLR              cm              1e+18                0               None False   True
      R_BLR_out                 BLR              cm              2e+18                0               None False   True
         L_Disk                Disk         erg / s              1e+45                0               None False  False
     R_inner_Sw                Disk      Sw. radii*                3.0                0               None False  False
       R_ext_Sw                Disk      Sw. radii*              500.0                0               None False  False
         T_Disk                Disk               K           100000.0                0               None False  False
       accr_eff                Disk                               0.08                0               None False  False
      disk_type                Disk                                 BB             None               None False   True
           M_BH                Disk          M_sun*       1000000000.0                0               None False  False
    -------------------------------------------------------------------------------------------------------------------


change Disk type
~~~~~~~~~~~~~~~~

the disk type can be set a mono temperature BB (as in the default case)
or as a more realistic multi temperature BB

.. code:: ipython3

    my_jet.set_par('disk_type',val='MultiBB')


now we set some parameter for the model

.. code:: ipython3

    my_jet.set_par('L_Disk',val=1E46)
    my_jet.set_par('gmax',val=5E4)
    my_jet.set_par('gmin',val=2.)
    my_jet.set_par('R_H',val=3E17)
    
    my_jet.set_par('p',val=1.5)
    my_jet.set_par('p_1',val=3.2)
    my_jet.set_par('R',val=3E15)
    my_jet.set_par('B',val=1.5)
    my_jet.set_par('z_cosm',val=0.6)
    my_jet.set_par('BulkFactor',val=20)
    my_jet.set_par('theta',val=1)
    my_jet.set_par('gamma_break',val=5E2)
    my_jet.set_N_from_nuLnu(nu_src=3E13,nuLnu_src=5E45)
    my_jet.set_IC_nu_size(100)

.. code:: ipython3

    my_jet.eval()
    p=my_jet.plot_model(frame='obs')
    p.rescale(y_min=-13.5,y_max=-9.5,x_min=9,x_max=27)



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_12_0.png


Dusty Torus
~~~~~~~~~~~

.. code:: ipython3

    my_jet.add_EC_component('DT')
    my_jet.show_model()


.. parsed-literal::

    
    -------------------------------------------------------------------------------------------------------------------
    jet model description
    -------------------------------------------------------------------------------------------------------------------
    name: BLR example  
    
    electron distribution:
     type: bkn  
     electron energy grid size:  1001
     gmin grid : 2.000000e+00
     gmax grid : 5.000000e+04
     normalization  True
     log-values  False
    
    radiative fields:
     seed photons grid size:  100
     IC emission grid size:  100
     source emissivity lower bound :  1.000000e-120
     spectral components:
       name:Sum, state: on
       name:Sync, state: self-abs
       name:SSC, state: on
       name:EC_BLR, state: on
       name:Disk, state: on
       name:EC_Disk, state: on
       name:DT, state: on
    external fields transformation method: blob
    
    SED info:
     nu grid size :200
     nu mix (Hz): 1.000000e+06
     nu max (Hz): 1.000000e+30
    
    flux plot lower bound   :  1.000000e-30
    
        name          par type           units             val         phys. bound. min  phys. bound. max   log  frozen
    ----------- ------------------- --------------- ------------------ ---------------- ------------------ ----- ------
              N    electron_density         1 / cm3  4174.081522033596                0               None False  False
           gmin  low-energy-cut-off lorentz-factor*                2.0                1       1000000000.0 False  False
           gmax high-energy-cut-off lorentz-factor*            50000.0                1 1000000000000000.0 False  False
              p   LE_spectral_slope                                1.5            -10.0               10.0 False  False
            p_1   HE_spectral_slope                                3.2            -10.0               10.0 False  False
    gamma_break    turn-over-energy lorentz-factor*              500.0                1       1000000000.0 False  False
              R         region_size              cm 3000000000000000.0           1000.0              1e+30 False  False
            R_H     region_position              cm              3e+17                0               None False   True
              B      magnetic_field               G                1.5                0               None False  False
          theta   jet-viewing-angle             deg                  1                0               None False  False
     BulkFactor     jet-bulk-factor Lorentz-factor*                 20              1.0               None False  False
         z_cosm            redshift                                0.6                0               None False  False
        tau_BLR                 BLR                                0.1                0                1.0 False  False
       R_BLR_in                 BLR              cm              1e+18                0               None False   True
      R_BLR_out                 BLR              cm              2e+18                0               None False   True
         L_Disk                Disk         erg / s              1e+46                0               None False  False
     R_inner_Sw                Disk      Sw. radii*                3.0                0               None False  False
       R_ext_Sw                Disk      Sw. radii*              500.0                0               None False  False
         T_Disk                Disk               K           100000.0                0               None False  False
       accr_eff                Disk                               0.08                0               None False  False
      disk_type                Disk                            MultiBB             None               None False   True
           M_BH                Disk          M_sun*       1000000000.0                0               None False  False
           T_DT                  DT               K              100.0                0               None False  False
           R_DT                  DT              cm              5e+18                0               None False  False
         tau_DT                  DT                                0.1                0                1.0 False  False
    -------------------------------------------------------------------------------------------------------------------


.. code:: ipython3

    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-13.5,y_max=-9.5,x_min=9,x_max=27)



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_15_0.png


.. code:: ipython3

    my_jet.add_EC_component('EC_DT')
    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-13.5,y_max=-9.5,x_min=9,x_max=27)



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_16_0.png


Changing the external field transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Default method, is the transformation of the external photon field from
the disk/BH frame to the relativistic blob

.. code:: ipython3

    my_jet.set_external_field_transf('blob')

Alternatively, in the case of istropric fields as the CMB or the BLR and
DT within the BLR radius, and DT radius, respectively, the it is
possible to transform the the electron distribution, moving the blob to
the disk/BH frame.

.. code:: ipython3

    my_jet.set_external_field_transf('disk')

External photon field energy density along the jet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def iso_field_transf(L,R,BulckFactor):
        beta=1.0 - 1/(BulckFactor*BulckFactor)
        return L/(4*np.pi*R*R*3E10)*BulckFactor*BulckFactor*(1+((beta**2)/3))
    
    def external_iso_behind_transf(L,R,BulckFactor):
        beta=1.0 - 1/(BulckFactor*BulckFactor)
        return L/((4*np.pi*R*R*3E10)*(BulckFactor*BulckFactor*(1+beta)**2))


EC seed photon fields, in the Disk rest frame

.. code:: ipython3

    %matplotlib inline
    fig = plt.figure(figsize=(8,6))
    ax=fig.subplots(1)
    N=50
    G=1
    R_range=np.logspace(13,25,N)
    y=np.zeros((8,N))
    my_jet.set_verbosity(0)
    my_jet.set_par('R_BLR_in',1E17)
    my_jet.set_par('R_BLR_out',1.1E17)
    for ID,R in enumerate(R_range):
        my_jet.set_par('R_H',val=R)
        my_jet.set_external_fields()
        my_jet.energetic_report(verbose=False)
        
        y[1,ID]=my_jet.energetic_dict['U_BLR_DRF']
        y[0,ID]=my_jet.energetic_dict['U_Disk_DRF']
        y[2,ID]=my_jet.energetic_dict['U_DT_DRF']
        
    y[4,:]=iso_field_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_DT.val,my_jet.parameters.R_DT.val,G)
    y[3,:]=iso_field_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_BLR.val,my_jet.parameters.R_BLR_in.val,G)
    y[5,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_BLR.val,R_range,G)
    y[6,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_DT.val,R_range,G)
    y[7,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative,R_range,G)
    
    ax.plot(np.log10(R_range),np.log10(y[0,:]),label='Disk')
    ax.plot(np.log10(R_range),np.log10(y[1,:]),'-',label='BLR')
    ax.plot(np.log10(R_range),np.log10(y[2,:]),label='DT')
    ax.plot(np.log10(R_range),np.log10(y[3,:]),'--',label='BLR uniform')
    ax.plot(np.log10(R_range),np.log10(y[4,:]),'--',label='DT uniform')
    ax.plot(np.log10(R_range),np.log10(y[5,:]),'--',label='BLR 1/R2')
    ax.plot(np.log10(R_range),np.log10(y[6,:]),'--',label='DT 1/R2')
    ax.plot(np.log10(R_range),np.log10(y[7,:]),'--',label='Disk 1/R2')
    ax.set_xlabel('log(R_H) cm')
    ax.set_ylabel('log(Uph) erg cm-3 s-1')
    
    ax.legend()





.. parsed-literal::

    <matplotlib.legend.Legend at 0x1826721390>




.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_25_1.png


.. code:: ipython3

    %matplotlib inline
    
    fig = plt.figure(figsize=(8,6))
    ax=fig.subplots(1)
    
    L_Disk=1E45
    N=50
    G=my_jet.parameters.BulkFactor.val
    R_range=np.logspace(15,22,N)
    y=np.zeros((8,N))
    my_jet.set_par('L_Disk',val=L_Disk)
    my_jet._blob.theta_n_int=100
    my_jet._blob.l_n_int=100
    my_jet._blob.theta_n_int=100
    my_jet._blob.l_n_int=100
    for ID,R in enumerate(R_range):
        my_jet.set_par('R_H',val=R)
        my_jet.set_par('R_BLR_in',1E17*(L_Disk/1E45)**.5)
        my_jet.set_par('R_BLR_out',1.1E17*(L_Disk/1E45)**.5)
        my_jet.set_par('R_DT',2.5E18*(L_Disk/1E45)**.5)
        my_jet.set_external_fields()
        my_jet.energetic_report(verbose=False)
        
        y[1,ID]=my_jet.energetic_dict['U_BLR']
        y[0,ID]=my_jet.energetic_dict['U_Disk']
        y[2,ID]=my_jet.energetic_dict['U_DT']
        
    
    
    y[4,:]=iso_field_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_DT.val,my_jet.parameters.R_DT.val,G)
    y[3,:]=iso_field_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_BLR.val,my_jet.parameters.R_BLR_in.val,G)
    y[5,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_BLR.val,R_range,G)
    y[6,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative*my_jet.parameters.tau_DT.val,R_range,G)
    y[7,:]=external_iso_behind_transf(my_jet._blob.L_Disk_radiative,R_range,G)
    
    ax.plot(np.log10(R_range),np.log10(y[0,:]),label='Disk')
    ax.plot(np.log10(R_range),np.log10(y[1,:]),'-',label='BLR')
    ax.plot(np.log10(R_range),np.log10(y[2,:]),'-',label='DT')
    ax.plot(np.log10(R_range),np.log10(y[3,:]),'--',label='BLR uniform')
    ax.plot(np.log10(R_range),np.log10(y[4,:]),'--',label='DT uniform')
    ax.plot(np.log10(R_range),np.log10(y[5,:]),'--',label='BLR 1/R2')
    ax.plot(np.log10(R_range),np.log10(y[6,:]),'--',label='DT 1/R2')
    ax.plot(np.log10(R_range),np.log10(y[7,:]),'--',label='Disk 1/R2')
    ax.axvline(np.log10( my_jet.parameters.R_DT.val ))
    ax.axvline(np.log10( my_jet.parameters.R_BLR_out.val))
    
    ax.set_xlabel('log(R_H) cm')
    ax.set_ylabel('log(Uph`) erg cm-3 s-1')
    
    ax.legend()





.. parsed-literal::

    <matplotlib.legend.Legend at 0x1825da1c10>




.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_26_1.png


IC against the CMB
~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    my_jet=Jet(name='test_equipartition',electron_distribution='lppl',beaming_expr='bulk_theta')
    my_jet.set_par('R',val=1E21)
    my_jet.set_par('z_cosm',val= 0.651)
    my_jet.set_par('B',val=2E-5)
    my_jet.set_par('gmin',val=50)
    my_jet.set_par('gamma0_log_parab',val=35.0E3)
    my_jet.set_par('gmax',val=30E5)
    my_jet.set_par('theta',val=12.0)
    my_jet.set_par('BulkFactor',val=3.5)
    my_jet.set_par('s',val=2.58)
    my_jet.set_par('r',val=0.42)
    my_jet.set_N_from_nuFnu(5E-15,1E12)
    my_jet.add_EC_component('EC_CMB')

We can now compare the different beaming pattern for the EC emission if
the CMB, and realize that the beaming pattern is different. This is very
important in the case of radio galaxies. The ``src`` transformation is
the one to use in the case of radio galaies or misaligned AGNs, and
gives a more accurate resut. Anyhow,be careful that this works onlyt for
isotropic external fields, suchs as the CMB, of the Disk and BLR seed
photons whitin the Dusty torus radius, and BLR radius, respectively

.. code:: ipython3

    from jetset.plot_sedfit import PlotSED
    p=PlotSED()
    
    my_jet.set_external_field_transf('blob')
    c= ['k', 'g', 'r', 'c'] 
    for ID,theta in enumerate(np.linspace(2,20,4)):
        my_jet.parameters.theta.val=theta
        my_jet.eval()
        my_jet.plot_model(plot_obj=p,comp='Sum',label='blob, theta=%2.2f'%theta,line_style='--',color=c[ID])
    
    my_jet.set_external_field_transf('disk')
    for ID,theta in enumerate(np.linspace(2,20,4)):
        my_jet.parameters.theta.val=theta
        my_jet.eval()
        my_jet.plot_model(plot_obj=p,comp='Sum',label='disk, theta=%2.2f'%theta,line_style='',color=c[ID])
    
    p.rescale(y_min=-17.5,y_max=-12.5,x_max=28)



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_30_0.png


Equipartition
-------------

It is also possible to set our jet at the equipartition, that is
achieved not using analytical approximation, but by numerically finding
the equipartition value over a grid. We have to provide the value of the
observed flux (``nuFnu_obs``) at a given observed frequency
(``nu_obs``), the minimum value of B (``B_min``), and the number of grid
points (``N_pts``)

.. code:: ipython3

    my_jet.parameters.theta.val=12
    B_min,b_grid,U_B,U_e=my_jet.set_B_eq(nuFnu_obs=5E-15,nu_obs=1E12,B_min=1E-9,N_pts=50,plot=True)
    my_jet.show_pars()
    
    my_jet.eval()
    p=my_jet.plot_model()
    p.rescale(y_min=-16.5,y_max=-13.5,x_max=28)


.. parsed-literal::

    B grid min  1e-09
    B grid max  1.0
    grid points 50



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_33_1.png


.. parsed-literal::

    setting B to  0.0001389495494373139
    setting N to  9.13927847193837e-06
          name             par type           units               val          phys. bound. min  phys. bound. max   log  frozen
    ---------------- ------------------- --------------- --------------------- ---------------- ------------------ ----- ------
                   N    electron_density         1 / cm3  9.13927847193837e-06              0.0               None False  False
                gmin  low-energy-cut-off lorentz-factor*                  50.0              1.0       1000000000.0 False  False
                gmax high-energy-cut-off lorentz-factor*             3000000.0              1.0 1000000000000000.0 False  False
                   s   LE_spectral_slope                                  2.58            -10.0               10.0 False  False
                   r  spectral_curvature                                  0.42            -15.0               15.0 False  False
    gamma0_log_parab    turn-over-energy lorentz-factor*               35000.0              1.0       1000000000.0 False  False
                   R         region_size              cm                 1e+21           1000.0              1e+30 False  False
                 R_H     region_position              cm                 1e+17              0.0               None False   True
                   B      magnetic_field               G 0.0001389495494373139              0.0               None False  False
               theta   jet-viewing-angle             deg                  12.0              0.0               None False  False
          BulkFactor     jet-bulk-factor Lorentz-factor*                   3.5              1.0               None False  False
              z_cosm            redshift                                 0.651              0.0               None False  False



.. image:: Jet_example_phys_EC_files/Jet_example_phys_EC_33_3.png


.. bibliography:: references.bib
