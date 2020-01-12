Numerical setup {#jet_numerical_guide}
===============

Changing the grid size for the electron distribution
----------------------------------------------------

``` {.sourceCode .ipython3}
from jetset.jet_model import Jet
my_jet=Jet(name='test',electron_distribution='lppl',)
my_jet.show_model()
```

::: {.parsed-literal}
### jet model description

name: test

electron distribution:

:   type: lppl electron energy grid size: 1001 gmin grid : 2.000000e+00
    gmax grid : 1.000000e+06 normalization True log-values False

radiative fields:

:   seed photons grid size: 100 IC emission grid size: 50 source
    emissivity lower bound : 1.000000e-120 spectral components:
    name:Sum, state: on name:Sync, state: self-abs name:SSC, state: on

external fields transformation method: blob

SED info:

:   nu grid size :200 nu mix (Hz): 1.000000e+06 nu max (Hz):
    1.000000e+30

flux plot lower bound : 1.000000e-30

> name par type units val phys. bound. min phys. bound. max log frozen

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False

:   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0 False
    False gmax high-energy-cut-off lorentz-factor\* 1000000.0 1.0
    1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0 10.0
    False False r spectral\_curvature 0.4 -15.0 15.0 False False

gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False

:   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False

R\_H region\_position cm 1e+17 0.0 None False True

:   B magnetic\_field G 0.1 0.0 None False False

beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False

:   z\_cosm redshift 0.1 0.0 None False False

------------------------------------------------------------------------
:::

It is possible to change the size of the grid for the electron
distributions. It is worth noting that at lower values of the grid size
the speed will increase, **but it is not recommended to go below 100**.

The actual value of the grid size is returned by the
`.Jet.gamma_grid_size`{.interpreted-text role="meth"}

``` {.sourceCode .ipython3}
print (my_jet.gamma_grid_size)
```

::: {.parsed-literal}
1001
:::

and this value can be changed using the method
`.Jet.set_gamma_grid_size`{.interpreted-text role="meth"}. In the
following we show the result for a grid of size=10, as anticipated the
final integration will be not satisfactory

``` {.sourceCode .ipython3}
my_jet.set_gamma_grid_size(10)
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_8_0.png)

``` {.sourceCode .ipython3}
my_jet.set_gamma_grid_size(100)
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_9_0.png)

``` {.sourceCode .ipython3}
my_jet.set_gamma_grid_size(1000)
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_10_0.png)

Changing the grid size for the seed photons
-------------------------------------------

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',)
my_jet.show_model()
```

::: {.parsed-literal}
### jet model description

name: test

electron distribution:

:   type: lppl electron energy grid size: 1001 gmin grid : 2.000000e+00
    gmax grid : 1.000000e+06 normalization True log-values False

radiative fields:

:   seed photons grid size: 100 IC emission grid size: 50 source
    emissivity lower bound : 1.000000e-120 spectral components:
    name:Sum, state: on name:Sync, state: self-abs name:SSC, state: on

external fields transformation method: blob

SED info:

:   nu grid size :200 nu mix (Hz): 1.000000e+06 nu max (Hz):
    1.000000e+30

flux plot lower bound : 1.000000e-30

> name par type units val phys. bound. min phys. bound. max log frozen

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False

:   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0 False
    False gmax high-energy-cut-off lorentz-factor\* 1000000.0 1.0
    1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0 10.0
    False False r spectral\_curvature 0.4 -15.0 15.0 False False

gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False

:   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False

R\_H region\_position cm 1e+17 0.0 None False True

:   B magnetic\_field G 0.1 0.0 None False False

beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False

:   z\_cosm redshift 0.1 0.0 None False False

------------------------------------------------------------------------
:::

we can get the current value of the seed photons grid size using
attribute `.Jet.nu_seed_size`{.interpreted-text role="meth"}

**in the current version there is lit of the size to 1000**

``` {.sourceCode .ipython3}
print (my_jet.nu_seed_size)
```

::: {.parsed-literal}
100
:::

and this value can be changed using the method
`.Jet.set_seed_nu_size`{.interpreted-text role="meth"}. In the following
we show the result for a grid of nu\_size=10

``` {.sourceCode .ipython3}
my_jet.nu_seed_size=10
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_17_0.png)

Changing the grid size for the IC process spectra
-------------------------------------------------

**in the current version there is a limit of the size to 1000**

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',)
my_jet.show_model()
```

::: {.parsed-literal}
### jet model description

name: test

electron distribution:

:   type: lppl electron energy grid size: 1001 gmin grid : 2.000000e+00
    gmax grid : 1.000000e+06 normalization True log-values False

radiative fields:

:   seed photons grid size: 100 IC emission grid size: 50 source
    emissivity lower bound : 1.000000e-120 spectral components:
    name:Sum, state: on name:Sync, state: self-abs name:SSC, state: on

external fields transformation method: blob

SED info:

:   nu grid size :200 nu mix (Hz): 1.000000e+06 nu max (Hz):
    1.000000e+30

flux plot lower bound : 1.000000e-30

> name par type units val phys. bound. min phys. bound. max log frozen

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False

:   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0 False
    False gmax high-energy-cut-off lorentz-factor\* 1000000.0 1.0
    1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0 10.0
    False False r spectral\_curvature 0.4 -15.0 15.0 False False

gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False

:   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False

R\_H region\_position cm 1e+17 0.0 None False True

:   B magnetic\_field G 0.1 0.0 None False False

beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False

:   z\_cosm redshift 0.1 0.0 None False False

------------------------------------------------------------------------
:::

``` {.sourceCode .ipython3}
print(my_jet.IC_nu_size)
```

::: {.parsed-literal}
50
:::

``` {.sourceCode .ipython3}
my_jet.IC_nu_size=20
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_22_0.png)

``` {.sourceCode .ipython3}
my_jet.IC_nu_size=100
my_jet.eval()
sed_plot=my_jet.plot_model()
sed_plot.rescale(x_min=8,y_min=-20,y_max=-12)
```

![image](Jet_example_num_files/Jet_example_num_23_0.png)
