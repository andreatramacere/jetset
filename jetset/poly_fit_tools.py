__author__ = "Andrea Tramacere"


import numpy as np

from jetset.poly_fit import do_log_Parab_FIT
from jetset.loglog_poly_model import LogLinear
from jetset.model_manager import FitModel
from jetset.minimizer import fit_SED
from jetset.data_loader import Data,ObsData


def get_SED_pl_fit(j, comp, fit_range):
    """
    power-law fit of a jet SED spectral component
    Parameters
    ----------
    j: jet object
    comp: spectral component name
    fit_range: [nu1,nu2] range (log Hz)

    Returns
    -------
    best-fit-par
    best-fit-err
    pl-func

    """

    c = j.get_spectral_component_by_name(comp)
    x = c.SED.nu.value
    y = c.SED.nuFnu.value

    data = Data(n_rows=x.shape[0])
    data.set_field('x', x)
    data.set_field('y', y)
    # data.set_field('dy',value=1E-15)
    data.set_meta_data('z', j.parameters.z_cosm.val)
    data.set_meta_data('restframe', 'obs')
    data.set_meta_data('data_scale', 'lin-lin')

    sed_data = ObsData(data_table=data)
    loglog_poly = LogLinear()
    loglog_pl = FitModel(name=None, loglog_poly=loglog_poly)
    mm, best_fit = fit_SED(loglog_pl,
                           sed_data,
                           10 ** fit_range[0],
                           10 ** fit_range[1],
                           loglog=True,
                           silent=True,
                           fitname='spectral-indices-best-fit',
                           minimizer='lsb')

    par = loglog_pl.get_par_by_name(loglog_poly, 'alpha')

    return par.best_fit_val, par.best_fit_err, loglog_pl


def get_SED_log_par_fit(x_p,y_p,j,comp,delta_p=[-1,1]):
    """
    log-parabolic fit of a  jet SED spectral component

    Parameters
    ----------
    x_p: starting value (log)
    y_p: starting value (log)
    j: jet object
    comp
    delta_p: +/- boundaries centered on x_p fit range (log of nu in Hz), i.e. fit rage= [log10(nu_p) + delta_p[0], nu_p + log10(delta_p[1]) ]

    Returns
    -------
    p [x_p, y_p, curvature]
    err [err x_p, err y_p, err curvature]

    """
    c=j.get_spectral_component_by_name(comp)
    msk = c.SED.nuFnu.value > 0
    x = np.log10(c.SED.nu.value[msk])
    y = np.log10(c.SED.nuFnu.value[msk])
    #x=np.log10(c.SED.nu.value)
    #y=np.log10(c.SED.nuFnu.value)

    p,err=do_log_Parab_FIT(x,y,x_p,y_p,-0.1,x_range=[x_p + delta_p[0], x_p + delta_p[1] ],dy=np.ones(x.size))
    p,err=do_log_Parab_FIT(x,y,p[0],p[1],p[2],x_range=[p[0] + delta_p[0], p[0] + delta_p[1] ],dy=np.ones(x.size))

    return p, err


def get_nu_p_S_delta_approx(my_jet,gp):
    """
    estimates nu peak of synchrotron delta-approx

    Parameters
    ----------
    my_jet: jet object
    gp: gamma-peak

    Returns
    -------
    nu_p

    """
    B=my_jet.parameters.B.val
    delta=my_jet.get_beaming()
    z=my_jet.parameters.z_cosm.val
    return np.log10(3.2E6*B*delta/(1+z))+2*gp


def get_n_gamma_log_par_fit(n, power=0, delta_p=[-0.5, 0.5]):
    """
    estimates the peak of n(gamma)*gamma^3 by log-par fit
    Parameters
    ----------
    n: emitters distribution  object
    power: 0 for n(gamma), 1 for n(gamma)*gamma, 2 for n(gamma)*gamma^2 etc...
    delta_p: +/- boundaries centered on x_p fit range (log), i.e. fit rage= [log10(gamma_p) + delta_p[0], log10(gamma_p) + log10(delta_p[1]) ]


    Returns
    -------
    p [x_p, y_p, curvature]
    err [err x_p, err y_p, err curvature]
    """

    n.update()

    if n.emitters_type == 'electrons' :
        x = np.log10(n.gamma_e)
        y = np.log10(n.n_gamma_e)
    elif n.emitters_type == 'protons':
        x = np.log10(n.gamma_p)
        y = np.log10(n.n_gamma_p)
    else:
        raise RuntimeError('emitters must be protons or electrons')

    y = y + power * x

    y_p = y.max()
    x_p = x[np.argmax(y)]
    p, err = do_log_Parab_FIT(x, y, x_p, y_p, -1, x_range=[x_p + delta_p[0], x_p + delta_p[1]], dy=np.ones(x.size))
    return p,err
