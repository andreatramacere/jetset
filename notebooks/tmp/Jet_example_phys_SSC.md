physical setup {#jet_physical_guide}
==============

In this section we describe how to build a model of jet able to
reproduce SSC/EC emission processes, using the `.Jet`{.interpreted-text
role="class"} class from the `.jet_model`{.interpreted-text role="mod"}
module. to This class through a flexible and intuitive interface allows
to access the C numerical code that provides an accurate and fast
computation of the synchrotron and inverse Compton fdsfsd processes.

basic setup
-----------

A jet instance can be built using the the `.Jet`{.interpreted-text
role="class"} class, istanciating the object in the following way:

``` {.sourceCode .ipython3}
from jetset.jet_model import Jet
my_jet=Jet(name='test',electron_distribution='lppl',)
```

This instruction will create:

:   -   a `Jet` object with `name` **test**,
    -   using as electron distribution the **lppl** model, that is a
        log-parabola with a low-energy power-law branch.
    -   using as working directory **test\_jet\_prod**

For a list of possible distribution you can run the command

``` {.sourceCode .ipython3}
Jet.available_electron_distributions()
```

::: {.parsed-literal}
lp: log-parabola pl: powerlaw lppl: log-parabola with low-energy
powerlaw branch lpep: log-parabola defined by peak energy plc: powerlaw
with cut-off bkn: broken powerlaw spitkov: spitkov lppl\_pile\_up:
log-parabola with low-energy powerlaw branch and pile-up bkn\_pile\_up:
broken powerlaw and pileup
:::

to view all the paramters:

``` {.sourceCode .ipython3}
my_jet.show_pars()
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False
>
> :   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0
>     False False gmax high-energy-cut-off lorentz-factor\* 1000000.0
>     1.0 1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0
>     10.0 False False r spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False
>
> :   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False
>
> R\_H region\_position cm 1e+17 0.0 None False True
>
> :   B magnetic\_field G 0.1 0.0 None False False
>
> beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False
>
> :   z\_cosm redshift 0.1 0.0 None False False
>
Each parameter has default values. All the parameters listed are handled
by `.ModelParameterArray`{.interpreted-text role="class"}, and each
parameter is an instance of the the `.JetParameter`{.interpreted-text
role="class"}. class. These parameters are also accessible as an astropy
table, with units:

``` {.sourceCode .ipython3}
my_jet.parameters.par_table
```

<i>Table length=11</i>
<table id="table47838080336" class="table-striped table-bordered table-condensed">
<thead><tr><th>name</th><th>par type</th><th>units</th><th>val</th><th>phys. bound. min</th><th>phys. bound. max</th><th>log</th><th>frozen</th></tr></thead>
<thead><tr><th>str16</th><th>str19</th><th>object</th><th>float64</th><th>float64</th><th>object</th><th>bool</th><th>bool</th></tr></thead>
<tr><td>N</td><td>electron_density</td><td>1 / cm3</td><td>100.0</td><td>0.0</td><td>None</td><td>False</td><td>False</td></tr>
<tr><td>gmin</td><td>low-energy-cut-off</td><td>lorentz-factor*</td><td>2.0</td><td>1.0</td><td>1000000000.0</td><td>False</td><td>False</td></tr>
<tr><td>gmax</td><td>high-energy-cut-off</td><td>lorentz-factor*</td><td>1000000.0</td><td>1.0</td><td>1000000000000000.0</td><td>False</td><td>False</td></tr>
<tr><td>s</td><td>LE_spectral_slope</td><td></td><td>2.0</td><td>-10.0</td><td>10.0</td><td>False</td><td>False</td></tr>
<tr><td>r</td><td>spectral_curvature</td><td></td><td>0.4</td><td>-15.0</td><td>15.0</td><td>False</td><td>False</td></tr>
<tr><td>gamma0_log_parab</td><td>turn-over-energy</td><td>lorentz-factor*</td><td>10000.0</td><td>1.0</td><td>1000000000.0</td><td>False</td><td>False</td></tr>
<tr><td>R</td><td>region_size</td><td>cm</td><td>5000000000000000.0</td><td>1000.0</td><td>1e+30</td><td>False</td><td>False</td></tr>
<tr><td>R_H</td><td>region_position</td><td>cm</td><td>1e+17</td><td>0.0</td><td>None</td><td>False</td><td>True</td></tr>
<tr><td>B</td><td>magnetic_field</td><td>G</td><td>0.1</td><td>0.0</td><td>None</td><td>False</td><td>False</td></tr>
<tr><td>beam_obj</td><td>beaming</td><td>Lorentz-factor*</td><td>10.0</td><td>0.0001</td><td>None</td><td>False</td><td>False</td></tr>
<tr><td>z_cosm</td><td>redshift</td><td></td><td>0.1</td><td>0.0</td><td>None</td><td>False</td><td>False</td></tr>
</table>
this means that you can easily convert the values in the table using the
units module of astropy.

::: {.warning}
::: {.admonition-title}
Warning
:::

Please note, that the table is built on the fly from the
`.ModelParameterArray`{.interpreted-text role="class"} and each
modification you do to this table will not be reflected on the actual
parameters values
:::

To get a full description of the model you can use the instruction

``` {.sourceCode .ipython3}
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

as you can notice, you can now access further information regarding the
model, such as numerical configuration of the grida. These parameters
will be discussed in the :ref:\`jet\_numerical\_guide\' section

If you want to use a comoslogy model different from the dafault one
please read the `cosmology`{.interpreted-text role="ref"} section.

::: {.warning}
::: {.admonition-title}
Warning
:::

Starting from version 1.1.0, the [R]{.title-ref} parameter as default is
linear and not logarithmic, please update your old scripts setting
[R]{.title-ref} with linear values.
:::

setting the parameters
----------------------

assume you want to change some of the parameters in your model, you can
use two methods:

1)  using the `.Jet.set_par()`{.interpreted-text role="class"} method

``` {.sourceCode .ipython3}
my_jet.set_par('B',val=0.2)
my_jet.set_par('gamma0_log_parab',val=5E3)
my_jet.set_par('gmin',val=1E2)
my_jet.set_par('gmax',val=1E8)
my_jet.set_par('R',val=1E15)
my_jet.set_par('N',val=1E3)
```

2)  accessing directly the parameter

``` {.sourceCode .ipython3}
my_jet.parameters.B.val=0.2
my_jet.parameters.r.val=0.4
```

investigating the electron distribution
---------------------------------------

``` {.sourceCode .ipython3}
my_jet.show_electron_distribution()
```

::: {.parsed-literal}

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--electron distribution:

:   type: lppl electron energy grid size: 1001 gmin grid : 2.000000e+00
    gmax grid : 1.000000e+06 normalization True log-values False

    > name par type units val phys. bound. min phys. bound. max log
    > frozen

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- B magnetic\_field G 0.2 0.0 None False False

:   N electron\_density 1 / cm3 1000.0 0.0 None False False R
    region\_size cm 1000000000000000.0 1000.0 1e+30 False False

R\_H region\_position cm 1e+17 0.0 None False True

:   beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False

gamma0\_log\_parab turn-over-energy lorentz-factor\* 5000.0 1.0 1000000000.0 False False

:   gmax high-energy-cut-off lorentz-factor\* 100000000.0 1.0
    1000000000000000.0 False False gmin low-energy-cut-off
    lorentz-factor\* 100.0 1.0 1000000000.0 False False r
    spectral\_curvature 0.4 -15.0 15.0 False False s LE\_spectral\_slope
    2.0 -10.0 10.0 False False

> z\_cosm redshift 0.1 0.0 None False False
:::

``` {.sourceCode .ipython3}
p=my_jet.electron_distribution.plot()
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_25_0.png)

``` {.sourceCode .ipython3}
p=my_jet.electron_distribution.plot3p()
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_26_0.png)

``` {.sourceCode .ipython3}
import numpy as np
p=None
for r in np.linspace(0.3,1,10):
    my_jet.parameters.r.val=r
    if p is None:
        p=my_jet.electron_distribution.plot3p()
    else:
        p=my_jet.electron_distribution.plot3p(p)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_27_0.png)

#### using log values for electron distribution parameters

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',electron_distribution_log_values=True)
my_jet.show_model()
```

::: {.parsed-literal}
### jet model description

name: test

electron distribution:

:   type: lppl electron energy grid size: 1001 gmin grid : 2.000000e+00
    gmax grid : 1.000000e+06 normalization True log-values True

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

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False

:   gmin low-energy-cut-off lorentz-factor\* 0.3010299956639812 0.0 9.0
    True False gmax high-energy-cut-off lorentz-factor\* 6.0 0.0 15.0
    True False s LE\_spectral\_slope 2.0 -10.0 10.0 False False r
    spectral\_curvature 0.4 -15.0 15.0 False False

gamma0\_log\_parab turn-over-energy lorentz-factor\* 4.0 0.0 9.0 True False

:   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False

R\_H region\_position cm 1e+17 0.0 None False True

:   B magnetic\_field G 0.1 0.0 None False False

beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False

:   z\_cosm redshift 0.1 0.0 None False False

------------------------------------------------------------------------
:::

evaluate and plot the model
---------------------------

At this point we can evaluate the emission for this jet model using the
instruction

``` {.sourceCode .ipython3}
my_jet.eval()
```

``` {.sourceCode .ipython3}
my_jet.show_pars()
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False
>
> :   gmin low-energy-cut-off lorentz-factor\* 0.3010299956639812 0.0
>     9.0 True False gmax high-energy-cut-off lorentz-factor\* 6.0 0.0
>     15.0 True False s LE\_spectral\_slope 2.0 -10.0 10.0 False False r
>     spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 4.0 0.0 9.0 True False
>
> :   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False
>
> R\_H region\_position cm 1e+17 0.0 None False True
>
> :   B magnetic\_field G 0.1 0.0 None False False
>
> beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False
>
> :   z\_cosm redshift 0.1 0.0 None False False
>
and plot the corresponding SED:

``` {.sourceCode .ipython3}
from jetset.plot_sedfit import PlotSED
my_plot=PlotSED()
my_jet.plot_model(plot_obj=my_plot)
my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_35_0.png)

alternatively, you can call the `plot_model` method without passing a
`Plot` object

``` {.sourceCode .ipython3}
my_plot=my_jet.plot_model()
my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_37_0.png)

If you want to have more points on the IC spectrum you can set the
numerical parameters for radiavite fields(see
:ref:\`jet\_numerical\_guide\' section for more details):

``` {.sourceCode .ipython3}
my_jet.set_IC_nu_size(100)
```

``` {.sourceCode .ipython3}
my_jet.eval()
my_plot=my_jet.plot_model()
my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_40_0.png)

you can access the same plot, but in the rest frame of the black hole,
or accretion disk, hence plotting the istropic luminosity, by simply
passing the `frame` kw to `src`

``` {.sourceCode .ipython3}
my_plot=my_jet.plot_model(frame='src')
my_plot.rescale(y_max=43,y_min=38,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_42_0.png)

the `my_plot` object returned will be built on the fly by the
`plot_model` method

if you wanto to have interacitve plot:

1)  in a jupyter notebook use:

``` {.sourceCode .no}
%matplotlib notebook
```

2) in jupyter lab:

:   ``` {.sourceCode .no}
    %matplotlib notebook
    ```

3)  in an ipython terminal

``` {.sourceCode .python}
from matplotlib import pylab as plt
plt.ion()
```

comparing models on the same plot
---------------------------------

to compare the same model after changing a parameter

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',)
my_jet.set_par('B',val=0.2)
my_jet.set_par('gamma0_log_parab',val=5E3)
my_jet.set_par('gmin',val=1E2)
my_jet.set_par('gmax',val=1E8)
my_jet.set_par('R',val=10**14.5)
my_jet.set_par('N',val=1E3)

my_jet.parameters.gamma0_log_parab.val=1E4
my_jet.eval()
my_plot=my_jet.plot_model(label='gamma0_log_parab=1E4',comp='Sum')
my_jet.set_par('gamma0_log_parab',val=1.0E5)
my_jet.eval()
my_plot=my_jet.plot_model(my_plot,label='gamma0_log_parab=1E5',comp='Sum')
my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_47_0.png)

saving a plot
-------------

to save the plot

``` {.sourceCode .ipython3}
my_plot.save('jet1.png')
```

saving and loading a model
--------------------------

::: {.warning}
::: {.admonition-title}
Warning
:::

starting from version 1.1.0 the saved model format has changed, if you
have models saved vith version\<1.1.0, plase update them the new models
by loading the old models with the
`.Jet.load_old_model`{.interpreted-text role="meth"} and then saving
them again.
:::

``` {.sourceCode .ipython3}
my_jet.save_model('test_model.dat')
```

``` {.sourceCode .ipython3}
my_jet_new=Jet.load_model('test_model.dat')
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- R region\_size cm 316227766016837.94 1000.0 1e+30 False False
>
> :   
>
>     R\_H region\_position cm 1e+17 0.0 None False True
>
>     :   B magnetic\_field G 0.2 0.0 None False False
>
> beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False
>
> :   
>
>     z\_cosm redshift 0.1 0.0 None False False
>
>     :   N electron\_density 1 / cm3 1000.0 0.0 None False False
>
>     gmin low-energy-cut-off lorentz-factor\* 100.0 1.0 1000000000.0 False False
>
>     :   
>
>         gmax high-energy-cut-off lorentz-factor\* 100000000.0 1.0 1000000000000000.0 False False
>
>         :   s LE\_spectral\_slope 2.0 -10.0 10.0 False False r
>             spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 100000.0 1.0
> 1000000000.0 False False

switching on/off the particle distribution normalization
--------------------------------------------------------

As default the electron distributions are normalized, i.e. are
multiplied by a constant `N_0`, in such a way that :

$\int_{\gamma_{min}}^{\gamma_{max}} n(\gamma) d\gamma =1$,

it means the the value [N]{.title-ref}, refers to the actual density of
emitters. If you want to chance this behavior, you can start looking at
the sate of `Norm_distr` flag with the following command

``` {.sourceCode .ipython3}
my_jet.Norm_distr
```

::: {.parsed-literal}
1
:::

and then you can switch off the normalization withe command

``` {.sourceCode .ipython3}
my_jet.switch_Norm_distr_OFF()
```

OR

``` {.sourceCode .ipython3}
my_jet.Norm_distr=0
```

or set back the normalization on with

``` {.sourceCode .ipython3}
my_jet.switch_Norm_distr_ON()
```

OR

``` {.sourceCode .ipython3}
my_jet.Norm_distr=1
```

setting the particle density from observed Fluxes or Luminosities
-----------------------------------------------------------------

It is possible to set the density of emitting particle starting from
some observed luminosity or flux (see the method
`.Jet.set_N_from_nuFnu`{.interpreted-text role="meth"},
meth:[.Jet.set\_N\_from\_nuLnu]{.title-ref})

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl')
```

this is the initial value of N

``` {.sourceCode .ipython3}
my_jet.parameters.N.val
```

::: {.parsed-literal}
100.0
:::

we now want to set the value of `N` in order that the observed
synchrotron flux at a given frequency matches a desired value. For
example, assume that we wish to set `N` in order that the synchrotron
flux at math:[10\^{15}]{.title-ref} Hz is exactly matching the desired
value of $10^{-=14}$ ergs cm-2 s-1. We can accomplish this by using the
`.Jet.get_par_by_name()`{.interpreted-text role="class"} as follows:

``` {.sourceCode .ipython3}
my_jet.set_N_from_nuFnu(nuFnu_obs=1E-14,nu_obs=1E15)
```

This is the updated value of `N`, obtained in order to match the given
flux at the given frequency

``` {.sourceCode .ipython3}
my_jet.get_par_by_name('N').val
```

::: {.parsed-literal}
271.77338679726074
:::

``` {.sourceCode .ipython3}
my_jet.parameters.show_pars()
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 271.77338679726074 0.0 None False False
>
> :   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0
>     False False gmax high-energy-cut-off lorentz-factor\* 1000000.0
>     1.0 1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0
>     10.0 False False r spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False
>
> :   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False
>
> R\_H region\_position cm 1e+17 0.0 None False True
>
> :   B magnetic\_field G 0.1 0.0 None False False
>
> beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False
>
> :   z\_cosm redshift 0.1 0.0 None False False
>
``` {.sourceCode .ipython3}
my_jet.eval()
my_plot=my_jet.plot_model(label='set N from F=1E-14')
my_plot.rescale(y_max=-13,y_min=-17.5,x_min=8)
```

![image](Jet_example_phys_SSC_files/Jet_example_phys_SSC_76_0.png)

as you can see, the synchrotron flux at $10^{15}$ Hz is exactly matching
the desired value of $10^{-14}$ ergs cm-2 s-1. Alternatively, the value
of N can be obtained using the rest-frame luminosity and frequency,
using the :class:\`.Jet.set\_N\_from\_nuLnu()

``` {.sourceCode .ipython3}
my_jet.set_N_from_nuLnu(nuLnu_src=1E43,nu_src=1E15)
```

where `L_0` is the source rest-frame istropic luminosity in erg/s at the
rest-frame frequency `nu_0` in Hz.

\#\# setting the beaming factor

It is possible to set the beaming factor according to the relativistic
BulkFactor and viewing angle, this can be done by setting the
`beaming_expr` kw in the Jet constructor, possible choices are

-   [delta]{.title-ref} to provide directly the beaming factor (default)
-   [bulk\_theta]{.title-ref} to provide the BulkFactor and the jet
    viewing angle

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')
```

``` {.sourceCode .ipython3}
my_jet.parameters.show_pars()
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False
>
> :   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0
>     False False gmax high-energy-cut-off lorentz-factor\* 1000000.0
>     1.0 1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0
>     10.0 False False r spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False
>
> :   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False
>
> R\_H region\_position cm 1e+17 0.0 None False True
>
> :   B magnetic\_field G 0.1 0.0 None False False
>
> theta jet-viewing-angle deg 0.1 0.0 None False False
>
> :   
>
>     BulkFactor jet-bulk-factor Lorentz-factor\* 10.0 1.0 None False False
>
>     :   z\_cosm redshift 0.1 0.0 None False False
>
the actual value of the beaming factor can be obtained using the
`.Jet.get_beaming`{.interpreted-text role="meth"}

``` {.sourceCode .ipython3}
my_jet.get_beaming()
```

::: {.parsed-literal}
19.943844732554165
:::

We can change the value of `theta` and get the updated value of the
beaming factor

``` {.sourceCode .ipython3}
my_jet.set_par('theta',val=10.)
```

``` {.sourceCode .ipython3}
my_jet.get_beaming()
```

::: {.parsed-literal}
4.968041140891955
:::

of course setting `beaming_expr=delta` we get the same beaming
expression as in the default case

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='delta')
```

``` {.sourceCode .ipython3}
my_jet.parameters.show_pars()
```

::: {.parsed-literal}
name par type units val phys. bound. min phys. bound. max log frozen
:::

> \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-- \-\-\-\-\-- N electron\_density 1 / cm3 100.0 0.0 None False False
>
> :   gmin low-energy-cut-off lorentz-factor\* 2.0 1.0 1000000000.0
>     False False gmax high-energy-cut-off lorentz-factor\* 1000000.0
>     1.0 1000000000000000.0 False False s LE\_spectral\_slope 2.0 -10.0
>     10.0 False False r spectral\_curvature 0.4 -15.0 15.0 False False
>
> gamma0\_log\_parab turn-over-energy lorentz-factor\* 10000.0 1.0 1000000000.0 False False
>
> :   R region\_size cm 5000000000000000.0 1000.0 1e+30 False False
>
> R\_H region\_position cm 1e+17 0.0 None False True
>
> :   B magnetic\_field G 0.1 0.0 None False False
>
> beam\_obj beaming Lorentz-factor\* 10.0 0.0001 None False False
>
> :   z\_cosm redshift 0.1 0.0 None False False
>
accessing individual spectral components
----------------------------------------

It is possible to access specific spectral components of our model

``` {.sourceCode .ipython3}
my_jet=Jet(name='test',electron_distribution='lppl',beaming_expr='bulk_theta')
my_jet.eval()
```

We can obtain this information anytime using the
`.Jet.list_spectral_components`{.interpreted-text role="meth"} method

``` {.sourceCode .ipython3}
my_jet.list_spectral_components()
```

::: {.parsed-literal}
Sum Sync SSC
:::

the on-screen message is telling us which components have been
evaluated.

and we cann access a specific component using the
`.Jet.get_spectral_component_by_name`{.interpreted-text role="meth"}
method

``` {.sourceCode .ipython3}
Sync=my_jet.get_spectral_component_by_name('Sync')
```

OR

``` {.sourceCode .ipython3}
Sync=my_jet.spectral_components.Sync
```

and from the `SED` object we can extract both the nu and nuFnu array

``` {.sourceCode .ipython3}
nu_sync=Sync.SED.nu
nuFnu_sync=Sync.SED.nuFnu
```

``` {.sourceCode .ipython3}
print (nuFnu_sync[::10])
```

::: {.parsed-literal}

\[1.00000000e-120 1.00000000e-120 1.08448642e-022 1.71738565e-018

:   4.07807919e-016 1.63686337e-015 6.48484725e-015 2.46700674e-014
    7.28812086e-014 1.24298363e-013 1.12162549e-013 1.42017250e-014
    4.14261886e-028 1.00000000e-120 1.00000000e-120 1.00000000e-120
    1.00000000e-120 1.00000000e-120 1.00000000e-120 1.00000000e-120\]
    erg / (cm2 s)
:::

or for the `src` rest frame (isotropic luminosity)

``` {.sourceCode .ipython3}
nu_sync_src=Sync.SED.nu_src
nuLnu_sync_src=Sync.SED.nuLnu_src
```

``` {.sourceCode .ipython3}
print (nuLnu_sync_src[::10])
```

::: {.parsed-literal}

\[2.70118406e-65 2.70118406e-65 2.92939742e+33 4.63897473e+37

:   1.10156425e+40 4.42146923e+40 1.75167660e+41 6.66383927e+41
    1.96865559e+42 3.35752756e+42 3.02971688e+42 3.83614730e+41
    1.11899760e+28 2.70118406e-65 2.70118406e-65 2.70118406e-65
    2.70118406e-65 2.70118406e-65 2.70118406e-65 2.70118406e-65\] erg /
    s
:::

Moreover, you can access the corresponding astropy table

``` {.sourceCode .ipython3}
my_jet.spectral_components.build_table(restframe='obs')
t_obs=my_jet.spectral_components.table
```

``` {.sourceCode .ipython3}
t_obs[::10]
```

<i>Table length=20</i>
<table id="table103693358800" class="table-striped table-bordered table-condensed">
<thead><tr><th>nu</th><th>Sum</th><th>Sync</th><th>SSC</th></tr></thead>
<thead><tr><th>Hz</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th><th>erg / (cm2 s)</th></tr></thead>
<thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>1000000.0</td><td>1e-120</td><td>1e-120</td><td>1e-120</td></tr>
<tr><td>15848931.924611142</td><td>1e-120</td><td>1e-120</td><td>1e-120</td></tr>
<tr><td>251188643.1509582</td><td>1.0844864302391386e-22</td><td>1.0844864159585182e-22</td><td>1.4280620238888498e-30</td></tr>
<tr><td>3981071705.5349693</td><td>1.7173856696253947e-18</td><td>1.7173856494785146e-18</td><td>2.0146879965841637e-26</td></tr>
<tr><td>63095734448.019424</td><td>4.078079552134786e-16</td><td>4.0780791901432327e-16</td><td>3.6199155354751893e-23</td></tr>
<tr><td>1000000000000.0</td><td>1.6368645001464381e-15</td><td>1.6368633684904028e-15</td><td>1.1316560354507247e-21</td></tr>
<tr><td>15848931924611.11</td><td>6.484856227306914e-15</td><td>6.484847252386628e-15</td><td>8.974920286013819e-21</td></tr>
<tr><td>251188643150958.22</td><td>2.467012104235855e-14</td><td>2.4670067379508708e-14</td><td>5.366284984149398e-20</td></tr>
<tr><td>3981071705534969.5</td><td>7.288148974605722e-14</td><td>7.288120857008097e-14</td><td>2.8117597624218914e-19</td></tr>
<tr><td>6.309573444801943e+16</td><td>1.2429968684418058e-13</td><td>1.2429836269824645e-13</td><td>1.324145934137875e-18</td></tr>
<tr><td>1e+18</td><td>1.1216821849060403e-13</td><td>1.1216254873133542e-13</td><td>5.669759268605211e-18</td></tr>
<tr><td>1.5848931924611109e+19</td><td>1.422429794065211e-14</td><td>1.420172495040777e-14</td><td>2.2572990244339145e-17</td></tr>
<tr><td>2.5118864315095718e+20</td><td>8.198273815038918e-17</td><td>4.142618855201174e-28</td><td>8.198273814997492e-17</td></tr>
<tr><td>3.9810717055349854e+21</td><td>2.6806229698492253e-16</td><td>1e-120</td><td>2.6806229698492253e-16</td></tr>
<tr><td>6.309573444801943e+22</td><td>7.79329160185085e-16</td><td>1e-120</td><td>7.79329160185085e-16</td></tr>
<tr><td>1e+24</td><td>1.876892626829062e-15</td><td>1e-120</td><td>1.876892626829062e-15</td></tr>
<tr><td>1.584893192461111e+25</td><td>2.722149689253548e-15</td><td>1e-120</td><td>2.722149689253548e-15</td></tr>
<tr><td>2.511886431509572e+26</td><td>9.2717312629558e-16</td><td>1e-120</td><td>9.2717312629558e-16</td></tr>
<tr><td>3.9810717055349856e+27</td><td>1e-120</td><td>1e-120</td><td>1e-120</td></tr>
<tr><td>6.309573444801943e+28</td><td>1e-120</td><td>1e-120</td><td>1e-120</td></tr>
</table>
and also in the `src` restframe

``` {.sourceCode .ipython3}
my_jet.spectral_components.build_table(restframe='src')
t_src=my_jet.spectral_components.table
```

``` {.sourceCode .ipython3}
t_src[::10]
```

<i>Table length=20</i>
<table id="table103693302224" class="table-striped table-bordered table-condensed">
<thead><tr><th>nu</th><th>Sum</th><th>Sync</th><th>SSC</th></tr></thead>
<thead><tr><th>Hz</th><th>erg / s</th><th>erg / s</th><th>erg / s</th></tr></thead>
<thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>1100000.0</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td></tr>
<tr><td>17433825.11707226</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td></tr>
<tr><td>276307507.4660541</td><td>2.929397454400055e+33</td><td>2.929397415825471e+33</td><td>3.8574583700258197e+25</td></tr>
<tr><td>4379178876.088467</td><td>4.638974788937108e+37</td><td>4.638974734516677e+37</td><td>5.442043094354434e+29</td></tr>
<tr><td>69405307892.82137</td><td>1.1015643465663552e+40</td><td>1.101564248785774e+40</td><td>9.77805812879182e+32</td></tr>
<tr><td>1100000000000.0</td><td>4.421472289763414e+40</td><td>4.421469232952174e+40</td><td>3.056811239929309e+34</td></tr>
<tr><td>17433825117072.223</td><td>1.751679024719035e+41</td><td>1.7516766004278765e+41</td><td>2.424291158119413e+35</td></tr>
<tr><td>276307507466054.06</td><td>6.663853762125038e+41</td><td>6.663839266801599e+41</td><td>1.449532343958061e+36</td></tr>
<tr><td>4379178876088467.0</td><td>1.9686631808560795e+42</td><td>1.9686555857754397e+42</td><td>7.595080639789025e+36</td></tr>
<tr><td>6.9405307892821384e+16</td><td>3.3575633227957894e+42</td><td>3.3575275551769376e+42</td><td>3.5767618852200225e+37</td></tr>
<tr><td>1.1000000000000001e+18</td><td>3.029870033860255e+42</td><td>3.0297168832268734e+42</td><td>1.5315063338183771e+38</td></tr>
<tr><td>1.743382511707222e+19</td><td>3.842244680626013e+41</td><td>3.836147300491401e+41</td><td>6.0973801346120284e+38</td></tr>
<tr><td>2.7630750746605293e+20</td><td>2.2145046516583798e+39</td><td>1.118997600209717e+28</td><td>2.21450465164719e+39</td></tr>
<tr><td>4.3791788760884844e+21</td><td>7.240856026525909e+39</td><td>2.7011840560827467e-65</td><td>7.240856026525909e+39</td></tr>
<tr><td>6.940530789282138e+22</td><td>2.1051115019323087e+40</td><td>2.7011840560827467e-65</td><td>2.1051115019323087e+40</td></tr>
<tr><td>1.1e+24</td><td>5.069832438569926e+40</td><td>2.7011840560827467e-65</td><td>5.069832438569926e+40</td></tr>
<tr><td>1.7433825117072222e+25</td><td>7.353027338882288e+40</td><td>2.7011840560827467e-65</td><td>7.353027338882288e+40</td></tr>
<tr><td>2.7630750746605295e+26</td><td>2.5044652659780157e+40</td><td>2.7011840560827467e-65</td><td>2.5044652659780157e+40</td></tr>
<tr><td>4.379178876088485e+27</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td></tr>
<tr><td>6.940530789282138e+28</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td><td>2.7011840560827467e-65</td></tr>
</table>
Of cousrse, since these colums have units, you can easily convert the
units of the Synchrotron luminostity form erg/s to GeV/s

``` {.sourceCode .ipython3}
t_src['Sync'][::10].to('GeV/s')
```

$$[1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.8283861 \times 10^{36},~2.8954203 \times 10^{40},~6.8754233 \times 10^{42},~2.759664 \times 10^{43},~1.0933105 \times 10^{44},~4.1592413 \times 10^{44},~1.2287382 \times 10^{45},~2.0956039 \times 10^{45},~1.8910005 \times 10^{45},~2.3943348 \times 10^{44},~6.9842337 \times 10^{30},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62},~1.6859465 \times 10^{-62}] \; \mathrm{\frac{GeV}{s}}$$

the table can be easily saved as an ascii file

``` {.sourceCode .ipython3}
t_src.write('test_SED.txt',format='ascii.ecsv',overwrite='True')
```

or in fits format

``` {.sourceCode .ipython3}
t_src.write('test_SED.fits',format='fits',overwrite='True')
```

::: {.parsed-literal}
WARNING: VerifyWarning: Keyword name \'restframe\' is greater than 8
characters or contains characters not allowed by the FITS standard; a
HIERARCH card will be created. \[astropy.io.fits.card\]
:::

Energetic report
----------------

It is possible to get an energetic report of the jet model (updated each
time that you eval the model). This report gives energy densities (`U_`)
(both in the blob end disk restframe), the luminosities of the emitted
components in the blob resftrame (`L_`), and the luminosity carried by
the jet (`jet_L`) for the radiative components, the electrons, the
magnetic fields, and for the cold protons in the jet.

``` {.sourceCode .ipython3}
my_jet.energetic_report()
```

::: {.parsed-literal}

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--jet eneregetic report:

:   name type units val

\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-- \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- U\_e Energy dens. blob rest. frame erg / cm3 0.0017404342430246782

:   U\_p Energy dens. blob rest. frame erg / cm3 0.015032764261 U\_B
    Energy dens. blob rest. frame erg / cm3 0.00039788735772973844

U\_Synch Energy dens. blob rest. frame erg / cm3 5.506769532122052e-05

:   

    U\_Synch\_DRF Energy dens. disk rest. frame erg / cm3 8.712292317747346

    :   

        U\_Disk Energy dens. blob rest. frame erg / cm3 0.0

        :   

            U\_BLR Energy dens. blob rest. frame erg / cm3 0.0

            :   U\_DT Energy dens. blob rest. frame erg / cm3 0.0

            U\_CMB Energy dens. blob rest. frame erg / cm3 0.0

    U\_Disk\_DRF Energy dens. disk rest. frame erg / cm3 0.0

    :   

        U\_BLR\_DRF Energy dens. disk rest. frame erg / cm3 0.0

        :   U\_DT\_DRF Energy dens. disk rest. frame erg / cm3 0.0

        U\_CMB\_DRF Energy dens. disk rest. frame erg / cm3 0.0
        L\_Sync\_rf Lum. blob rest. frme. erg / s 1.728764352592126e+38
        L\_SSC\_rf Lum. blob rest. frme. erg / s 3.82887909757934e+36

L\_EC\_Disk\_rf Lum. blob rest. frme. erg / s 0.0

:   

    L\_EC\_BLR\_rf Lum. blob rest. frme. erg / s 0.0

    :   L\_EC\_DT\_rf Lum. blob rest. frme. erg / s 0.0

    L\_EC\_CMB\_rf Lum. blob rest. frme. erg / s 0.0

    :   L\_PP\_rf Lum. blob rest. frme. erg / s 0.0

    jet\_L\_Sync jet Lum. erg / s 4.3219108814803147e+39

    :   jet\_L\_SSC jet Lum. erg / s 9.572197743948349e+37

jet\_L\_EC\_Disk jet Lum. erg / s 0.0

:   

    jet\_L\_EC\_BLR jet Lum. erg / s 0.0

    :   jet\_L\_EC\_DT jet Lum. erg / s 0.0

    jet\_L\_EC\_CMB jet Lum. erg / s 0.0

    :   jet\_L\_PP jet Lum. erg / s 0.0

    jet\_L\_rad jet Lum. erg / s 4.417632858919798e+39

    :   jet\_L\_kin jet Lum. erg / s 4.043042849486075e+42 jet\_L\_tot
        jet Lum. erg / s 4.047460482344995e+42 jet\_L\_e jet Lum. erg /
        s 4.097964612089291e+41 jet\_L\_B jet Lum. erg / s
        9.368514312500004e+40 jet\_L\_p jet Lum. erg / s
        3.539561245152146e+42

------------------------------------------------------------------------
:::

If you want to evaluate the energetic report in non verbose mode:

``` {.sourceCode .ipython3}
my_jet.energetic_report(verbose=False)
```

``` {.sourceCode .ipython3}
my_jet.energetic_dict
```

::: {.parsed-literal}

{\'U\_e\': 0.0017404342430246782,

:   \'U\_p\': 0.015032764261, \'U\_B\': 0.00039788735772973844,
    \'U\_Synch\': 5.506769532122052e-05, \'U\_Synch\_DRF\':
    8.712292317747346, \'U\_Disk\': 0.0, \'U\_BLR\': 0.0, \'U\_DT\':
    0.0, \'U\_CMB\': 0.0, \'U\_Disk\_DRF\': 0.0, \'U\_BLR\_DRF\': 0.0,
    \'U\_DT\_DRF\': 0.0, \'U\_CMB\_DRF\': 0.0, \'L\_Sync\_rf\':
    1.728764352592126e+38, \'L\_SSC\_rf\': 3.82887909757934e+36,
    \'L\_EC\_Disk\_rf\': 0.0, \'L\_EC\_BLR\_rf\': 0.0,
    \'L\_EC\_DT\_rf\': 0.0, \'L\_EC\_CMB\_rf\': 0.0, \'L\_PP\_rf\': 0.0,
    \'jet\_L\_Sync\': 4.3219108814803147e+39, \'jet\_L\_SSC\':
    9.572197743948349e+37, \'jet\_L\_EC\_Disk\': 0.0,
    \'jet\_L\_EC\_BLR\': 0.0, \'jet\_L\_EC\_DT\': 0.0,
    \'jet\_L\_EC\_CMB\': 0.0, \'jet\_L\_PP\': 0.0, \'jet\_L\_rad\':
    4.417632858919798e+39, \'jet\_L\_kin\': 4.043042849486075e+42,
    \'jet\_L\_tot\': 4.047460482344995e+42, \'jet\_L\_e\':
    4.097964612089291e+41, \'jet\_L\_B\': 9.368514312500004e+40,
    \'jet\_L\_p\': 3.539561245152146e+42}
:::

``` {.sourceCode .ipython3}
my_jet.energetic_report_table
```

<i>Table length=33</i>
<table id="table103711276688" class="table-striped table-bordered table-condensed">
<thead><tr><th>name</th><th>type</th><th>units</th><th>val</th></tr></thead>
<thead><tr><th>str13</th><th>str29</th><th>object</th><th>float64</th></tr></thead>
<tr><td>U_e</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.0017404342430246782</td></tr>
<tr><td>U_p</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.015032764261</td></tr>
<tr><td>U_B</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.00039788735772973844</td></tr>
<tr><td>U_Synch</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>5.506769532122052e-05</td></tr>
<tr><td>U_Synch_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>8.712292317747346</td></tr>
<tr><td>U_Disk</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.0</td></tr>
<tr><td>U_BLR</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.0</td></tr>
<tr><td>U_DT</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.0</td></tr>
<tr><td>U_CMB</td><td>Energy dens. blob rest. frame</td><td>erg / cm3</td><td>0.0</td></tr>
<tr><td>U_Disk_DRF</td><td>Energy dens. disk rest. frame</td><td>erg / cm3</td><td>0.0</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>jet_L_EC_BLR</td><td>jet Lum.</td><td>erg / s</td><td>0.0</td></tr>
<tr><td>jet_L_EC_DT</td><td>jet Lum.</td><td>erg / s</td><td>0.0</td></tr>
<tr><td>jet_L_EC_CMB</td><td>jet Lum.</td><td>erg / s</td><td>0.0</td></tr>
<tr><td>jet_L_PP</td><td>jet Lum.</td><td>erg / s</td><td>0.0</td></tr>
<tr><td>jet_L_rad</td><td>jet Lum.</td><td>erg / s</td><td>4.417632858919798e+39</td></tr>
<tr><td>jet_L_kin</td><td>jet Lum.</td><td>erg / s</td><td>4.043042849486075e+42</td></tr>
<tr><td>jet_L_tot</td><td>jet Lum.</td><td>erg / s</td><td>4.047460482344995e+42</td></tr>
<tr><td>jet_L_e</td><td>jet Lum.</td><td>erg / s</td><td>4.097964612089291e+41</td></tr>
<tr><td>jet_L_B</td><td>jet Lum.</td><td>erg / s</td><td>9.368514312500004e+40</td></tr>
<tr><td>jet_L_p</td><td>jet Lum.</td><td>erg / s</td><td>3.539561245152146e+42</td></tr>
</table>
