

__author__ = "Andrea Tramacere"


import matplotlib as mpl


try:
    from matplotlib import  pyplot as plt
except:
    try:
        from matplotlib import pylab as plt

    except:
        try:
           import  pylab as plt
        except:
            raise RuntimeError('Unable to import pylab/pyplot from matplotlib')



from matplotlib import gridspec
import numpy as np
import  os
from astropy.constants import m_e,m_p,c
import matplotlib.ticker as ticker
import warnings

from collections import namedtuple

from .output import section_separator,WorkPlace

from .utils import *

__all__=['PlotSED','BasePlot','PlotPdistr','PlotSpecComp','PlotSeedPhotons','PlotSpectralMultipl','PlotTempEvDiagram','PlotTempEvEmitters']

def y_ev_transf(x):
    return x / 2.417E14

def y_ev_transf_inv(x):
    return x * 2.417E14



def set_mpl():
    mpl.rcParams['figure.figsize'] = [12.0, 8.0]
    mpl.rcParams['figure.dpi'] = 100
    mpl.rcParams['savefig.dpi'] = 100

    mpl.rcParams['font.size'] = '14'
    mpl.rcParams['legend.fontsize'] = 'medium'
    mpl.rcParams['figure.titlesize'] = 'medium'



def _rescale( x_min=None, x_max=None, y_min=None, y_max=None):
        warnings.warn('`The rescale method has been removed and has been replaced by the setlim method')
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!The rescale method as been replaced by the setlim method            !")
        print("!please notice that now jetset uses log axis rather than loglog plots!")
        print("!so, the correct way to use it is rescale(x_min=8)->setlim(x_min=1E8)!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

class  PlotSED (object):
    def __init__(self,
                 sed_data=None,
                 model=None,
                 interactive=False,
                 plot_workplace=None,
                 title='Plot',
                 frame='obs',
                 density=False,
                 dpi=100,
                 figsize=(12,8),
                 use_grid=True):

        check_frame(frame)

        self.frame=frame

        self.axis_kw=['x_min','x_max','y_min','y_max']
        self.interactive=interactive

        plot_workplace=plot_workplace
        self.lines_data_list=[]
        self.lines_model_list=[]
        self.lines_res_list = []

        if self.interactive is True:
            plt.ion()
            print ('running PyLab in interactive mode')

        if plot_workplace is None:
            plot_workplace=WorkPlace()
            self.out_dir=plot_workplace.out_dir
            self.flag=plot_workplace.flag

        else:
            self.out_dir=plot_workplace.out_dir
            self.flag=plot_workplace.flag


            self.title="%s_%s"%(title,self.flag)

        if figsize is None:
            figsize=(10,8)

        self.fig=plt.figure(figsize=figsize,dpi=dpi)

        self.gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])

        self.sedplot= self.fig.add_subplot(self.gs[0])
        self._add_res_plot()

        self.set_plot_axis_labels(density=density)

        #if autoscale==True:
        self.sedplot.set_autoscalex_on(True)
        self.sedplot.set_autoscaley_on(True)
        self.sedplot.set_autoscale_on(True)
        self.counter=0
        if use_grid is True:
            self.sedplot.grid(use_grid,alpha=0.5)

        self.sedplot.set_xlim(1E6, 1E30)
        if frame == 'obs':
            self.sedplot.set_ylim(1E-20, 1E-8)
        elif frame == 'src':
            self.sedplot.set_ylim(1E38, 1E55)
        elif frame == 'blob':
            self.sedplot.set_ylim(1E34, 1E51)
        else:
            unexpected_behaviour()

        self.sedplot.set_xscale("log", nonpositive='clip')
        self.sedplot.set_yscale("log", nonpositive='clip')
   
        self.secaxy = self.sedplot.secondary_xaxis('top', functions=(y_ev_transf, y_ev_transf_inv))
        self.secaxy.set_xlabel('E (eV)')

        self.resplot.set_ybound(-2,2)
        
        try:
            if hasattr(self.fig.canvas.manager,'toolbar'):
                self.fig.canvas.manager.toolbar.update()
        except:
            pass

        if sed_data is not None :
            self.add_data_plot(sed_data,density=density)

        if model is not  None:
            self.add_model_plot(model)
        self.counter_res=0

        self.add_residual_plot = self.add_model_residual_plot

    def _check_frame(self,frame):
        if frame is None:
            frame=self.frame

        elif frame != self.frame:
            raise RuntimeError('you have to use the same restframe of the PlotSED class:',self.frame )

        return frame

    def _add_res_plot(self):
        self.resplot = self.fig.add_subplot(self.gs[1], sharex=self.sedplot)

        self.lx_res = '$ \\nu $  (Hz)'
        self.ly_res = 'res'

        self.resplot.set_ylabel(self.ly_res)
        self.resplot.set_xlabel(self.lx_res)
        #self,resplot.set_xscale("log", nonpositive='clip')
        self.add_res_zeroline()

    def clean_residuals_lines(self):
        for i in range(len(self.lines_res_list)):
            self.del_residuals_line(0)

    def clean_data_lines(self):

        for i in range(len(self.lines_data_list)):
            self.del_data_line(0)

    def clean_model_lines(self):
        for i in range(len(self.lines_model_list)):
            self.del_model_line(0)


    def list_lines(self):
        if self.lines_data_list==[] and self.lines_model_list==[]:
            pass
        else:

            for ID,plot_line in enumerate(self.lines_data_list):
                print('data',ID, plot_line.get_label())

            for ID,plot_line in enumerate(self.lines_model_list):
                print ('model',ID,  plot_line.get_label())

    def del_data_line(self,line_ID):
        if self.lines_data_list==[]:
            print  ("no lines to delete ")
        else:
            print ("removing line: ",self.lines_data_list[line_ID])
            line = self.lines_data_list[line_ID]

            for item in line:
                # This removes lines
                if np.shape(item) == ():
                    item.remove()
                else:
                    # This removes containers for data with errorbars
                    for item1 in item:
                        item1.remove()

            del self.lines_data_list[line_ID]
            #self.update_legend()
            self.update_plot()

    def del_model_line(self,line_ID):

        if self.lines_model_list==[]:
            #print  "no lines to delete "
            pass
        else:

            line=self.lines_model_list[line_ID]
            line.remove()


            del self.lines_model_list[line_ID]

            self.update_plot()
            #self.update_legend()

    def del_residuals_line(self, line_ID):
        if self.lines_res_list == []:
            # print  "no lines to delete "
            pass
        else:

            line = self.lines_res_list[line_ID]
            line.remove()

            del self.lines_res_list[line_ID]

            self.update_plot()
            #self.update_legend()

    def set_plot_axis_labels(self, density=False):
        self.lx = '$ \\nu $  (Hz)'

        if self.frame == 'src' or self.frame == 'blob':

            if density is False:
                self.ly = '$ \\nu L_{\\nu} $   (erg  s$^{-1})$'
            else:
                self.ly = '$   L_{\\nu} $   (erg  s$^{-1}$ Hz$^{-1})$'

        elif self.frame == 'obs':
            if density is False:
                self.ly = '$ \\nu F_{\\nu} $   (erg cm$^{-2}$  s$^{-1})$'
            else:
                    self.ly = '$   F{\\nu} $   (erg cm$^{-2}$  s$^{-1}$ Hz$^{-1})$'


        else:
            unexpected_behaviour()

        self.sedplot.set_ylabel(self.ly)
        self.sedplot.set_xlabel(self.lx)



    def add_res_zeroline(self):
        #y0 = np.zeros(2)
        #x0 = [0,30]
        self.resplot.axhline(0, ls='--', color='black')
        self.update_plot()

    
    def rescale(self, x_min=None, x_max=None, y_min=None, y_max=None):
        _rescale(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max)

    def setlim(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.sedplot.set_xlim(x_min, x_max)
        self.sedplot.set_ylim(y_min, y_max)
    
    def setlim_res(self,x_min=None,x_max=None,y_min=None,y_max=None):
        self.resplot.set_xlim(x_min,x_max)
        self.resplot.set_ylim(y_min,y_max)
        self.update_plot()

    def update_plot(self):
        self.fig.canvas.draw()
       
        y_s = []
        x_min = []
        x_max = []
        y_min = None
        y_max = None
        if len(self.sedplot.lines)>0:

            for l in self.sedplot.lines:
                if len(l.get_ydata())>0:
                    y_s.append(np.max(l.get_ydata()))
            if len(y_s) > 0:
                y_min = min(y_s)/1000
                y_max = max(y_s)*10
            else:
                self.sedplot.autoscale(axis='y')
            if y_min is not None and y_max is not None:
                self.sedplot.set_ylim(y_min, y_max)
                for l in self.sedplot.lines:
                    x=np.array(l.get_xdata())[np.array(l.get_ydata()) >= y_min]
                    if len(x)>0:
                        x_min.append(np.min(x))
                        x_max.append(np.max(x))
                if len(x_min)>0  and  len(x_max)>0:
                    self.sedplot.set_xlim(min(x_min)/10, max(x_max)*10)
        else:
            self.sedplot.relim()
            self.sedplot.autoscale(axis='y')
            self.sedplot.autoscale(axis='x')
    
        self.update_legend()
        self.fig.tight_layout()

    def update_legend(self,label=None):

        _handles=[]

        if self.lines_data_list!=[] and self.lines_data_list is not None:
            _handles.extend(self.lines_data_list)

        if self.lines_model_list!=[] and self.lines_model_list is not None:
            _handles.extend(self.lines_model_list)

        for h in _handles[:]:
            if h._label is  None:
                _handles.remove(h)
            elif h._label.startswith('_line'):
                _handles.remove(h)
            else:
                 pass

        self.sedplot.legend(handles=_handles,loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1, prop={'size':10})




    def add_model_plot(self, model, label=None, color=None, line_style=None, flim=None,auto_label=True,fit_range=None,density=False, update=True, lw=1.0 ,frame=None):

        frame=self._check_frame(frame=frame)

        if hasattr(model,'get_model_points'):
            try:
                x, y = model.get_model_points(log_log=False, frame = self.frame)
            except Exception as e:
                raise RuntimeError('for model',model.name, "problem with get_model_points()",e)
        else:
            try:
                x, y = model.SED.get_model_points(log_log=False, frame = self.frame)
            except Exception as e:
                raise RuntimeError('for model',model.name, "problem with SED.get_model_points()",e)

        if density is True:
            y=y/x
        if line_style is None:
            line_style = '-'

        if label is None and auto_label is True:
            if model.name is not None:
                label = model.name
            else:
                label = 'line %d' % self.counter

        if flim is not None:
            msk=y>flim
            x=x[msk]
            y=y[msk]
        else:
            pass

        if fit_range is not None:
            msk1 = x < fit_range[1]
            msk2 = x > fit_range[0]

            x = x[msk1 * msk2]
            y = y[msk1 * msk2]

        line, = self.sedplot.plot(x, y, line_style, label=label,color=color,lw=lw)


        self.lines_model_list.append(line)

        if update is True:
            #self.update_legend()
            self.update_plot()

        self.counter += 1

    def plot_tempev_model(self,
                          temp_ev,
                          region,
                          comp='Sum',
                          frame=None,
                          t1=None,
                          t2=None,
                          time_slice=None,
                          time_slice_bin=None,
                          time=None,
                          time_bin=None,
                          density=False,
                          use_cached=False,
                          sed_data=None,
                          average=False):

        frame=self._check_frame(frame)


        if (time_slice is not None and time is not None):
            raise RuntimeError('you can to pass either the N-th time slice "time_slice", or the blob time in seconds "time" ')

        if t1 is None or t1 < region.time_sampled_emitters.time_blob[0]:
            t1 = region.time_sampled_emitters.time_blob[0]

        if t2 is None or t2 > region.time_sampled_emitters.time_blob[-1]:
            t2 = region.time_sampled_emitters.time_blob[-1]

        if time_slice is None:
            _time_slice = 0
        else:
            _time_slice = time_slice

        _time_slice_bin = time_slice_bin
        if time_slice_bin is None and time_slice is None:
            _time_slice_bin = 1


        #if time_slice is not None or time_bin is None:

        if time is not None and time_bin is not None:
            t_array = np.arange(t1, t2, time_bin)
            time_id_array=None
        elif time is not None and time_bin is None:
            t_array = np.array([time])
            time_id_array = None
        else:
            t_array, time_id_array = region.time_sampled_emitters._get_time_samples(time_slice=_time_slice,
                                                                                    time_slice_bin=_time_slice_bin)
            time_id_array = time_id_array[t_array <= t2]
            time_id_array = time_id_array[t_array >= t1]
            t_array = t_array[t_array <= t2]
            t_array = t_array[t_array >= t1]


        g = plt.cm.Greens(np.linspace(0.5, 1, t_array.size))
        r = plt.cm.Reds(np.linspace(0.5, 1, t_array.size))
        b = plt.cm.Blues(np.linspace(0.5, 1, t_array.size))
        for ID, t in enumerate(t_array):
            if time is not None:
                s = region.get_SED(comp, frame=frame, time=t, use_cached=use_cached, time_bin=time_bin,average=average)
            else:
                s = region.get_SED(comp, frame=frame, time_slice=time_id_array[ID], use_cached=use_cached,
                                   time_slice_bin=time_slice_bin,average=average)

            label = None
            ls = '-'
            color = r[ID]
            if temp_ev.custom_q_jnj_profile[temp_ev._get_time_slice_T_array(t)] > 0:
                color = g[ID]
                ls = '-'
                lw = 0.2
            if temp_ev.custom_acc_profile[temp_ev._get_time_slice_T_array(t)] > 0:
                color = b[ID]
                ls = '-'
                lw = 0.2

            if ID == 0:
                lw = 2
                ls = '--'
                label = 'start, t=%2.2e (s)' % t
                color = 'green'
            if ID == t_array.size - 1:
                lw = 2
                ls = '--'
                color = 'purple'
                label = 'stop, t=%2.2e (s)' % t

            self.add_model_plot(model=s, label=label, line_style=ls, color=color, update=False, lw=lw,
                                    auto_label=False,density=density)

        if sed_data is not None:
            self.add_data_plot(sed_data)

        self.update_plot()
        return


    def add_data_plot(self,sed_data,label=None,color=None,frame=None,fmt='o',ms=4,mew=0.5,fit_range=None, density = False):

        frame = self._check_frame(frame)
        try:
            x,y,dx,dy,=sed_data.get_data_points(log_log=False,frame=self.frame, density=density)
        except Exception as e:
            raise RuntimeError("!!! ERROR failed to get data points from", sed_data,e)


        if dx is None:
            dx=np.zeros(len(sed_data.data['nu_data']))


        if dy is None:
            dy=np.zeros(len(sed_data.data['nu_data']))

        UL = sed_data.data['UL']

        if label is None:
            if sed_data.obj_name is not None  :
                label=sed_data.obj_name
            else:
                label='line %d'%self.counter

        if fit_range is not None:
            msk1 = x < fit_range[1]
            msk2 = x > fit_range[0]

            x = x[msk1 * msk2]
            y = y[msk1 * msk2]
            dx= dx[msk1 * msk2]
            dy = dy[msk1 * msk2]
            UL=UL[msk1 * msk2]

        line = self.sedplot.errorbar(x, y, xerr=dx, yerr=dy, fmt=fmt
                                     , uplims=UL,label=label,ms=ms,mew=mew,color=color)

        self.lines_data_list.append(line)

        self.counter+=1
        #self.update_legend()
        self.update_plot()
    


    def add_xy_plot(self,x,y,label=None,color=None,line_style=None,autoscale=False):

        if line_style is None:
            line_style='-'


        if label is None:
            label='line %d'%self.counter

        line, = self.sedplot.plot(x, y, line_style,label=label)

        self.lines_model_list.append(line)

        self.counter+=1

        #self.update_legend()
        self.update_plot()



    def add_model_residual_plot(self, model, data, label=None, color=None, filter_UL=True, fit_range=None):
        if data is not None:
            x,y = model.get_residuals(log_log=False,data=data,filter_UL=filter_UL)
            self.add_xy_residual_plot(x=x, y=y, fit_range=fit_range, color=color)
        else:
            pass


    def add_xy_residual_plot(self, x, y, fit_range=None, color=None):
        if self.counter_res == 0:
            self.add_res_zeroline()
        if fit_range is not None:
            msk1 = x < fit_range[1]
            msk2 = x > fit_range[0]

            x = x[msk1 * msk2]
            y = y[msk1 * msk2]

        line = self.resplot.errorbar(x, y, yerr=np.ones(x.size), fmt='+', color=color)
        self.lines_res_list.append(line)
        self.counter_res += 1
        self.update_plot()


    def add_text(self,lines):
        self.PLT.focus(0,0)
        x_min, x_max = self.sedplot.get_xlim()
        y_min, y_max = self.sedplot.get_ylim()
        t=''
        for line in lines:
            t+='%s \\n'%line.strip()
        self.PLT.text(t,font=10,charsize=0.6,x=x_min-1.5,y=y_min-2.85)
        self.PLT.redraw()


    def save(self,filename=None):
        if filename is None:
            wd=self.out_dir
            filename = 'jetset_fig.png'

        else:
            wd=''

        outname = os.path.join(wd,filename)
        self.fig.savefig(outname)

    def show(self):
        self.fig.show()





class BasePlot(object):

    def __init__(self,figsize=(8,6),dpi=100):
        self.fig, self.ax = plt.subplots(figsize=figsize,dpi=dpi)

    def rescale(self, x_min=None, x_max=None, y_min=None, y_max=None):
        _rescale(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max)

    def setlim(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

    def update_plot(self):
        self.fig.canvas.draw()
        self.ax.relim()
        self.ax.autoscale(axis='y')
        self.ax.legend()
        self.fig.tight_layout()



class PlotSpectralMultipl(BasePlot):
    def __init__(self):
        super(PlotSpectralMultipl, self).__init__()

        secax = self.ax.secondary_xaxis('top', functions=(y_ev_transf, y_ev_transf_inv))
        secax.set_xlabel('E (eV)')


    def plot(self,nu,y,y_label,y_min=None,y_max=None,label=None,line_style=None,color=None):

        self.ax.plot(np.log10(nu), np.log10(y),label=label,ls=line_style,color=color)
        self.ax.set_xlabel(r'$ \nu $  (Hz)')
        self.ax.set_ylabel(y_label)
        self.ax.set_ylim(y_min, y_max)
        self.ax.legend()
        self.update_plot()




class  PlotPdistr (BasePlot):

    def __init__(self,figsize=(8,6),dpi=100,injection=False,loglog=True):
        super(PlotPdistr, self).__init__(figsize=figsize,dpi=dpi)
        self.loglog=loglog
        self.injection = injection

    def _set_variable(self,gamma,n_gamma,particle,energy_unit,pow=None):

        energy_plot=False
        if energy_unit == 'gamma':
            energy_name = '\gamma'
            energy_units=''
        else:
            energy_name='E'
            energy_units= '%s'%energy_unit
            energy_plot=True

        if  energy_plot is False:
            x=gamma
            y=n_gamma

        else:

            if particle=='electrons':
                x = gamma*(m_e*c*c).to(energy_unit).value
                y = n_gamma * 1.0/(m_e*c*c).to(energy_unit).value
            elif particle=='protons':
                x = gamma * (m_p * c * c).to(energy_unit).value
                y = n_gamma * 1.0 / (m_p * c * c).to(energy_unit).value
            else:
                raise  RuntimeError('particle ',particle, 'not implemented')

        m = y > 0
        x=np.copy(x)
        y=np.copy(y)

        if pow is not None:
            y[m] = y[m]* np.power( x[m], pow)

        if self.loglog is True:
            x[m] = np.log10( x[m])
            y[m] = np.log10(y[m])

        return x[m], y[m], energy_name,energy_units



    def _set_xy_label(self,energy_name,energy_units,pow):
        if energy_units != '':
            _e = '(%s)' % energy_units
        else:
            _e = ''

        if self.loglog is True:
            self.ax.set_xlabel(r'log($%s$)  %s' % (energy_name, _e))
        else:
            self.ax.set_xlabel(r'$%s$  %s' % (energy_name, _e))

        if energy_units != '':
            _e = '%s^{-1}' % energy_units
        else:
            _e = ''

        n_str = 'n($%s$)'%energy_name
        if pow is not None and pow!=0:
            n_str = 'n($%s$) $%s^{%d}$' % (energy_name,energy_name,pow)
            if energy_units != '':
                i=-1+pow
                if i==1:
                    _e='%s' %(energy_units)
                else:
                    _e='%s^{%d}' %(energy_units,i)
        if self.injection is False:

            if self.loglog is True:
                self.ax.set_ylabel(r'log(%s)   ($cm^{-3} %s$) ' % (n_str, _e))
            else:
                self.ax.set_ylabel(r'%s   ($cm^{-3} %s$) ' % (n_str, _e))
        else:
            if self.loglog is True:
                self.ax.set_ylabel(r'log(Q$_{inj}$($%s$))   ($cm^{-3} s^{-1} %s$)' % (energy_name,_e))
            else:
                self.ax.set_ylabel(r'Q$_{inj}$($%s$)   ($cm^{-3} s^{-1} %s$)' % (energy_name,_e))

    def _plot(self,x,y,c=None,lw=None,ls=None,label=None):
        if self.loglog is True:
            self.ax.plot(x, y, c=c, lw=lw, label=label,ls=ls)
        else:
            self.ax.loglog(x, y,c=c, lw=lw, label=label,ls=ls)

    def plot_distr(self,gamma,n_gamma,y_min=None,y_max=None,x_min=None,x_max=None,particle='electrons',energy_unit='gamma',label=None):

        x,y,energy_name,energy_units=self._set_variable(gamma,n_gamma,particle,energy_unit)

        if label is None:
            label=particle
        self._plot(x,y,label=label)
        self._set_xy_label(energy_name,energy_units,pow=None)
        self.update_plot()
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlim(x_min, x_max)


    def plot_distr2p(self, gamma, n_gamma, y_min=None, y_max=None, x_min=None, x_max=None,particle='electrons',energy_unit='gamma',label=None):
        if label is None:
            label=particle

        x, y, energy_name, energy_units = self._set_variable(gamma, n_gamma, particle, energy_unit,pow=2)
        self._plot(x,y,label=label)
        self._set_xy_label(energy_name, energy_units,pow=2)
        self.update_plot()
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlim(x_min, x_max)

    def plot_distr3p(self,gamma,n_gamma,y_min=None,y_max=None,x_min=None,x_max=None,particle='electrons',energy_unit='gamma', label=None):
        if label is None:
            label = particle

        x, y, energy_name, energy_units = self._set_variable(gamma, n_gamma, particle, energy_unit, pow=3)
        self._plot(x,y,label=label)
        self._set_xy_label(energy_name, energy_units,pow=3)
        self.update_plot()
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlim(x_min, x_max)


    def update_plot(self):
        self.fig.canvas.draw()
        self.ax.relim()
        self.ax.autoscale(axis='y')
        self.ax.autoscale(axis='x')
        self.ax.legend()
        self.fig.tight_layout()







class  PlotTempEvEmitters (PlotPdistr):

    def __init__(self,figsize=(8,6),dpi=100,loglog=True):
        super(PlotTempEvEmitters, self).__init__(figsize=figsize,dpi=dpi,loglog=loglog,)


    def _plot_distr(self,temp_ev,
                        region,particle='electrons',
                        energy_unit='gamma',
                        pow=None,
                        plot_Q_inj=True,
                        t1=None, 
                        t2=None,
                        #time_slice_bin=None,
                        #time_slice=None):
                    ):

        if t1 is None:
            t1=region.time_sampled_emitters.time_blob[0]

        if t2 is None:
            t2=region.time_sampled_emitters.time_blob[-1]
        
       # if time_slice_bin is None and time_slice is None:
       #     _time_slice_bin = 1

        t_array = np.linspace(t1, t2, region.time_sampled_emitters.time_blob.size,)

        ls = '-'
        lw = 0.2

        g = plt.cm.Greens(np.linspace(0.5, 1, t_array.size))
        r = plt.cm.Reds(np.linspace(0.5, 1, t_array.size))
        b = plt.cm.Blues(np.linspace(0.5, 1, t_array.size))

        n = region.time_sampled_emitters.n_gamma
        for ID, t in enumerate(t_array):
            label = None
            ls = '-'
            color = r[ID]

            if temp_ev.custom_q_jnj_profile[temp_ev._get_time_slice_T_array(t)] > 0:
                color = g[ID]
                ls = '-'
                lw = 0.2

            if temp_ev.custom_acc_profile[temp_ev._get_time_slice_T_array(t)] > 0:
                color = b[ID]
                ls = '-'
                lw = 0.2

            if ID == 0:
                lw = 2
                ls = '--'
                label = 'start, t=%2.2e (s)' % t
                color = 'green'
            if ID == t_array.size - 1:
                lw = 2
                ls = '--'
                color = 'purple'
                label = 'stop, t=%2.2e (s)' % t

            x, y, energy_name, energy_units = self._set_variable(region.time_sampled_emitters.gamma, n[ID], particle, energy_unit, pow=pow)
            self._plot(x,y,c=color,lw=lw,label=label,ls=ls)
        #x, y, energy_name, energy_units = self._set_variable(region.time_sampled_emitters.gamma, n[0], particle, #energy_unit, pow=pow)
        #self._plot(x, y, c='black', lw=2,label='Start sample')

        #x, y, energy_name, energy_units = self._set_variable(region.time_sampled_emitters.gamma, n[-1], particle, #energy_unit, pow=pow)
        #self._plot(x, y, c='blue', lw=2,label='Stop sample')
        self._set_xy_label(energy_name, energy_units,pow=pow)
        #TODO move to plot inj if iny is used in region
        if temp_ev.Q_inj is not None and (region._region_type == 'acc' or temp_ev._only_radiation is True):
            y = temp_ev.Q_inj.n_gamma_e * temp_ev.delta_t
            x = temp_ev.Q_inj.gamma_e
            if plot_Q_inj is True:
                if pow is not None:
                    y=y*np.power(x,pow)

                self._plot(x,y, c='red', lw=1, label='$Q_{inj}$ delta t')
        #print('==> d')

        self.ax.legend()

    def plot_distr(self, temp_ev, region='acc', energy_unit='gamma',plot_Q_inj=True,pow=None):
        self._plot_distr(temp_ev,region=region, particle='electrons',energy_unit=energy_unit,pow=pow,plot_Q_inj=plot_Q_inj)

    def plot_distr2p(self, temp_ev, region='acc', energy_unit='gamma',plot_Q_inj=True):
        self._plot_distr(temp_ev, region=region, particle='electrons',energy_unit=energy_unit,pow=2,plot_Q_inj=plot_Q_inj)

    def plot_distr3p(self, temp_ev, region='acc', energy_unit='gamma',plot_Q_inj=True):
        self._plot_distr(temp_ev, region=region, particle='electrons',energy_unit=energy_unit,pow=3,plot_Q_inj=plot_Q_inj)


class  PlotTempEvDiagram (object):

    def __init__(self,figsize=(8,6),dpi=100,expanding_region=False):
        if expanding_region is True:
            n_rows =4
        else:
            n_rows = 4
        self.fig, self.axs = plt.subplots(n_rows, 1, figsize=figsize, dpi=dpi, sharex=False)

    def plot(self,
             T_array,
             inj_profile,
             acc_profile,
             R_exp,
             B_exp,
             R_H_exp):


        self.axs[0].plot(T_array, acc_profile, label='Acc. start/stop', c='g')
        self.axs[0].set_ylim(0, 1.5)

        self.axs[1].plot(T_array, inj_profile, label='Inj. profile', c='b')
        self.axs[1].set_ylim(0,None)
        self.axs[1].set_xlabel('Time in blob frame (s)')
        self.axs[1].sharex = self.axs[0]
        #if expanding_region is True:
        self.axs[2].plot(np.log10(R_H_exp), np.log10(R_exp), label='Rad. region size', c='orange')
        self.axs[2].set_ylabel('log(R_rad) (cm)')
        self.axs[3].plot(np.log10(R_H_exp), np.log10(B_exp), label='B Rad. region', c='green')
        self.axs[3].set_xlabel('log(BH distance) in observer frame (cm)')
        self.axs[3].set_ylabel('log(B) (G)')
        self.axs[3].sharex = self.axs[2]

        for ax in self.axs:
            ax.legend()


        #self.fig.subplots_adjust(hspace=0)
        self.fig.tight_layout()

    def rescale(self, x_min=None, x_max=None, y_min=None, y_max=None):
        _rescale(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max)

    def setlim(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

class  PlotSpecComp (BasePlot):

    def __init__(self):
        super(PlotSpecComp, self).__init__()

    def plot(self,nu,nuFnu,y_min=None,y_max=None):

        self.ax.plot(np.log10(nu), np.log10(nuFnu))
        self.ax.set_xlabel(r'log($ \nu $)  (Hz)')
        self.ax.set_ylabel(r'log($ \nu F_{\nu} $ )  (erg cm$^{-2}$  s$^{-1}$)')
        self.ax.set_ylim(y_min, y_max)
        self.update_plot()



class  PlotSeedPhotons (BasePlot):

    def __init__(self):
        super(PlotSeedPhotons, self).__init__()

    def plot(self,nu,nuFnu,y_min=None,y_max=None):

        self.ax.plot(np.log10(nu), np.log10(nuFnu))
        self.ax.set_xlabel(r'log($ \nu $)  (Hz)')
        self.ax.set_ylabel(r'log(n )  (photons cm$^{-3}$  Hz$^{-1}$  ster$^{-1}$)')
        self.ax.set_ylim(y_min, y_max)
        self.update_plot()
