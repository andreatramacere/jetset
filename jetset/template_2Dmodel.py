
__author__ = "Andrea Tramacere"



from scipy import interpolate

import numpy as np

from astropy.units import Unit as u
from astropy.units import spectral
from astropy.table import Table

import os

from  .plot_sedfit import PlotSED,PlotSpectralMultipl

from .model_parameters import ModelParameter
from .base_model import  Model, MultiplicativeModel

__all__=['TemplateTable2D', 'EBLAbsorptionTemplate']



class TemplateTable2D(Model):
    """
    Class to handle spectral templates
    """
    def __init__(self,
                 x_values,
                 y_values,
                 z_values,
                 log_input_grid=False,
                 log_log_interp=False,
                 x_in_units=None,
                 y_in_units=None,
                 x_out_units='Hz',
                 y_out_units='erg cm-2 s-1',
                 z_in_units=None,
                 z_out_units=None,
                 nu_size=100,
                 name='TableModel2D',
                 zero=1E-100):
        """
        """
        super(TemplateTable2D, self).__init__(name=name)


        self._zero=zero

        self._log_log_interp = log_log_interp

        if log_input_grid is False:

            self.x_values = x_values
            self.y_values = y_values
            self.z_values = z_values
        else:
            self.z_values = np.power(10., z_values)
            self.x_values = np.power(10., x_values)
            self.y_values = np.power(10., y_values)

        if x_in_units is not None:
            self.x_values *= u(x_in_units).to(x_out_units, equivalencies=spectral())

        if hasattr(self.x_values,'value'):
            _x = self.x_values.value
        else:
            _x = self.x_values

        if y_in_units is not None:
            self.y_values *= u(y_in_units).to(y_out_units, equivalencies=spectral())

        if hasattr(self.y_values,'value'):
            _y = self.y_values.value
        else:
            _y = self.y_values

        if z_in_units is not None:
            self.z_values *= u(z_in_units).to(z_out_units, equivalencies=spectral())

        if hasattr(self.z_values, 'value'):
            _z = self.z_values.value
        else:
            _z = self.z_values

        if log_log_interp is True:
            self._scale='log-log'
            _x[_x<self._zero]=self._zero
            self._x_grid = np.log10(_x)
            self._y_grid = np.log10(_y)
            self._z_grid = np.log10(_z)

        else:
            self._scale='lin-lin'
            self._x_grid =  _x
            self._y_grid =  _y
            self._z_grid =  _z


        self.interp_func = interpolate.RectBivariateSpline(self._x_grid, self._y_grid, self._z_grid)

        self.nu_size = nu_size
        
        self.name = name
        
        self.model_type = 'table2D'

    def plot_model(self,plot_obj=None,clean=False,label=None,sed_data=None,color=None, density=False,frame='obs'):
        plot_obj=self._set_up_plot(plot_obj,sed_data,frame,density)

        if clean==True:
            plot_obj.clean_model_lines()

        if label is None:
            label=self.name

        plot_obj.add_model_plot(self.SED, line_style='-', label=label, flim=self.flux_plot_lim,color=color, density=density,frame=frame)

        return plot_obj

    def _func(self,x,y):

        if np.shape(x)==():
                x=np.array([x])

        if self._log_log_interp is True:
            x=np.log10(x)
            y=np.log10(y)
             
        model = self.interp_func(x,y)
        #
        model[model<self._zero]=self._zero

        if self._log_log_interp is True:
            return np.power(10., model)
        else:
            pass

        return model

    
    def eval(self,fill_SED=True,x=None,y=None,get_model=False,loglog=False):
        """
        Evaluates the Template for the current parameters values
        """    

        if x is None:
            x = np.copy(self.x_values)

        if loglog is True:
            x = np.power(10,x)

        if y is None:
            y = np.copy(self.y_values)

        if loglog is False:
            y = np.power(10,y)
        
        

        
        model=self._func(x,y)



        if get_model==True:
            if loglog==False:
            
                return model
            else:
                
                return np.log10(model)

        else:
            return None


class EBLAbsorptionTemplate(TemplateTable2D,MultiplicativeModel):
    """

    """
    def __init__(self,
                 redshift_array,
                 energy_array,
                 tau_values,
                 applied_model=None,
                 z=1.0,
                 template_name='tau_ebl_template',
                 nu_size=100):
        """

        Parameters
        ----------
        redshift_array
        energy_array
        tau_values
        applied_model
        z
        template_name
        nu_size
        """
        super(EBLAbsorptionTemplate, self).__init__(
            x_values=redshift_array,
            y_values=energy_array,
            z_values=tau_values,
            log_input_grid=False,
            log_log_interp=False,
            nu_size=nu_size,
            name=template_name
        )



        if z is None and applied_model is not None:
            self.parameters.add_par(applied_model.get_par_by_type('redshift'))
        elif applied_model is None and z is not None:
            self.parameters.add_par(ModelParameter(name='z_cosm',
                                                   par_type='redshift',
                                                   val=z,
                                                   val_min=0,
                                                   val_max=None,
                                                   units='',
                                                   frozen=True,
                                                   log=False))
        elif  z is not None and applied_model is not None:
            raise RuntimeError('either you provide a redshift value (z) or model with redshift parameter')

    @property
    def z(self):
        if hasattr(self,'_z'):
            if isinstance(self._z,ModelParameter):
                return self._z.val
            else:
                return self._z
        else:
            return  None

    @z.setter
    def z(self,z):
        self._z = z

    def apply_to(self,model):
        p_m = model.parameters.get_par_by_type('redshift')
        _p= self.get_par_by_name('redshift')
        _p.val=p_m.val

    def _check(self):
        pass

    @classmethod
    def from_name(cls,template_name,applied_model=None,z=1.0,nu_size=100):
        """

        Parameters
        ----------
        template_name
        applied_model
        z
        nu_size

        Returns
        -------

        """

        _Templates_dir = os.path.dirname(__file__) + '/ebl_data'

        _allowed_templates = ['Finke_2010', 'Dominguez_2010', 'Franceschini_2008']

        if template_name not in _allowed_templates:
            raise ValueError('template EBL model', template_name, 'not in allowdr',_allowed_templates)

        _template_name_dict = {}

        _template_name_dict['Finke_2010'] = 'tau_finke_2010.fits'
        _template_name_dict['Dominguez_2010'] = 'tau_dominguez_2010.fits'
        _template_name_dict['Franceschini_2008'] = 'tau_franceschini_2008.dat'

        file_path = os.path.join(_Templates_dir, _template_name_dict[template_name])

        if file_path.endswith('fits'):
            data = Table.read(file_path,format='fits')
        elif file_path.endswith('dat'):
            data = Table.read(file_path, format='ascii.ecsv')
        else:
            data = Table.read(file_path)

        y_values = data['energies']
        y_values = np.log10(y_values.to('Hz', equivalencies=spectral()).value)

        try:
            x_values = np.array(data.meta['REDSHIFT'], dtype=np.float)
        except:
            x_values = np.array(data.meta['redshift'], dtype=np.float)

        cn = [name for name in data.colnames if name != 'energies']
        z_values = np.array([data[n].data for n in cn])

        return cls(redshift_array=x_values,
                   energy_array=y_values,
                   tau_values=z_values,
                   template_name=template_name,
                   applied_model=applied_model,
                   z=z,
                   nu_size=nu_size)

    def eval(self, fill_SED=True, nu=None, get_model=False, loglog=False):
        """

        Parameters
        ----------
        fill_SED
        nu
        get_model
        loglog
        z

        Returns
        -------

        """

        out_model=None

        if nu is None:
            log_nu = np.copy(self.y_values)

        else:
            if loglog is True:
                log_nu = nu
            else:
                log_nu = np.log10(nu)


        z=self.parameters.get_par_by_name('z_cosm').val

        model = np.exp(-self._func(z,log_nu).T).flatten()

        if get_model == True:
            if loglog == False:
                out_model=model
            else:
                out_model=  np.log10(model)

        if fill_SED is True:
            if loglog is False:
                _nu=np.power(10.,log_nu)
            else:
                _nu=np.copy(log_nu)

            self.nu=_nu
            self.tau=model

        return out_model

    def plot_model(self, plot_obj=None,  label=None, line_style='-',color=None,frame='obs'):
        if plot_obj is None:
            plot_obj = PlotSpectralMultipl()

        if label is None:
            label = '%s'%self.name
            label += ' ,z=%5.5f'%self.parameters.get_par_by_name('z_cosm').val

        plot_obj.plot(nu=self.nu,y=self.tau, y_label=r'$ log10(\exp^{- \tau}) $', line_style=line_style, label=label, color=color)


        return plot_obj
