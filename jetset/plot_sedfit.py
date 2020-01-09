from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

NOPYLAB=True

try:
    
    import  pylab as plt
   
    from matplotlib import pylab as pp
    from matplotlib import  gridspec

except:

    NOPYLAB=True

    #print "pylab not found on this system"
    #print "install package, and/or update pythonpath"



from collections import namedtuple

from .output import section_separator,WorkPlace

from .utils import *

import numpy as np


__all__=['PlotSED']


class  PlotSED (object):
    
    
    def __init__(self,sed_data=None,
                 model=None,
                 #x_min=6.0,
                 #x_max=3.0,
                 #y_min=-20.0,
                 #y_max=-9.0,
                 interactive=False,
                 plot_workplace=None,
                 title='Plot',
                 frame='obs',
                 figsize=None):
                 #autoscale=True):

        check_frame(frame)

        self.frame=frame

        self.axis_kw=['x_min','x_max','y_min','y_max']
        self.interactive=interactive

        plot_workplace=plot_workplace
        #self.line_tuple=namedtuple('line',['label','ref'])
        self.lines_data_list=[]
        self.lines_model_list=[]
        self.lines_res_list = []

        #self.legend=[]
     
        

        if self.interactive==True:
       
            
            plt.ion()
            print ('running PyLab in interactive mode')
        

        #--------------------------------------------------------------

        
        #set workplace
        if plot_workplace is None:
            plot_workplace=WorkPlace()
            self.out_dir=plot_workplace.out_dir
            self.flag=plot_workplace.flag
     
        else:
            self.out_dir=plot_workplace.out_dir
            self.flag=plot_workplace.flag
        
        
            self.title="%s_%s"%(title,self.flag)
        
        
        #Build sedplot    
      
        if figsize is None:
            figsize=(10,8)

        self.fig=plt.figure(figsize=figsize)
        


        self.gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
            
        
        self.sedplot= self.fig.add_subplot(self.gs[0])
        self._add_res_plot()
        
        self.set_plot_axis_labels()
        
        #if autoscale==True:
        self.sedplot.set_autoscalex_on(True)
        self.sedplot.set_autoscaley_on(True)
        self.sedplot.set_autoscale_on(True)
        self.counter=0



        self.sedplot.grid(True,alpha=0.5)
       


        self.sedplot.set_xlim(6, 30)
        if frame == 'obs':
           self.sedplot.set_ylim(-20, -8)

        elif frame == 'src':
            self.sedplot.set_ylim(38, 55)
        else:
            unexpetced_behaviour()


        self.resplot.set_ybound(-2,2)
        try:
            if hasattr(self.fig.canvas.manager,'toolbar'):
                self.fig.canvas.manager.toolbar.update()
        except:
            pass

        if sed_data is not None :
            self.add_data_plot(sed_data)

        if model is not  None:
            self.add_model_plot(model)


            
        self.counter_res=0


    def _add_res_plot(self):




        self.resplot = self.fig.add_subplot(self.gs[1], sharex=self.sedplot)

        self.lx_res = 'log($ \\nu $)  (Hz)'
        self.ly_res = 'res'

        self.resplot.set_ylabel(self.ly_res)
        self.resplot.set_xlabel(self.lx_res)



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

                # self.sedplot.lines.remove(self.lines_list[line_ID].ref[0])

            #self.legend.remove(self.lines_data_list[line_ID].label)

            del self.lines_data_list[line_ID]

            self.update_legend()
            self.update_plot()

            #self.fig.canvas.draw()
    
    
    def del_model_line(self,line_ID):
        
        if self.lines_model_list==[]:
            #print  "no lines to delete "
            pass
        else:

            line=self.lines_model_list[line_ID]
            line.remove()


            del self.lines_model_list[line_ID]


            self.update_plot()
            self.update_legend()

    def del_residuals_line(self, line_ID):

        if self.lines_res_list == []:
            # print  "no lines to delete "
            pass
        else:

            line = self.lines_res_list[line_ID]
            line.remove()

            del self.lines_res_list[line_ID]

            self.update_plot()
            self.update_legend()

    def set_plot_axis_labels(self):
        self.lx = 'log($ \\nu $)  (Hz)'

        if self.frame == 'src':

            self.ly = 'log($ \\nu L_{\\nu} $ )  (erg  s$^{-1}$)'
            
        elif self.frame == 'obs':

                self.ly = 'log($ \\nu F_{\\nu} $ )  (erg cm$^{-2}$  s$^{-1}$)'
        else:
            unexpetced_behaviour()

        self.sedplot.set_ylabel(self.ly)
        self.sedplot.set_xlabel(self.lx)
        
        #self.sedplot.set_xlim(self.x_min,self.x_max)
        #self.sedplot.set_ylim(self.y_min,self.y_max)
        
    
    def add_res_zeroline(self):


        y0 = np.zeros(2)
        x0 = [0,30]

        self.resplot.plot(x0,y0,'--',color='black')
        self.update_plot()

        
        
     
     
     
     
    def rescale(self,x_min=None,x_max=None,y_min=None,y_max=None):


            
        self.sedplot.set_xlim(x_min,x_max)
        self.sedplot.set_ylim(y_min,y_max)


    
    #def autoscale(self):

        
        #self.sedplot.autoscale_view(tight=True)
    #    for l in self.sedplot.lines:

    #        x_min,x_max=self.sedplot.get_xlim()

    #        y_min,y_max=self.sedplot.get_ylim()
        
    #    self.sedplot.set_xticks(np.arange(int(x_min)-2,int(x_max)+2,1.0))
        
    #    self.sedplot.set_xlim(x_min-1,x_max+1)
        
    #    self.sedplot.set_ylim(y_min-1,y_max+1)
        
        
        
      #  if self.resplot is not None  :
            
            #self.resplot.autoscale_view(tight=True)
            
            #self.x_min_res=self.x_min-1
            #self.x_max_res=self.x_max+1

        #    y_min_res,y_max_res=self.resplot.get_ylim()
        #    x_min_res,x_max_res=self.resplot.get_xlim()

            
        #    self.resplot.set_xticks(np.arange(int(x_min_res)-2,int(x_max_res)+2,1.0))
            
        #    self.resplot.set_xlim(x_min_res,x_max_res)
        #    self.resplot.set_ylim(y_min_res,y_max_res)
            
            
            
        #self.update_plot()
        
        #self.sedplot.set_autoscale_on(False)
   
   
    
    def rescale_res(self,x_min=None,x_max=None,y_min=None,y_max=None):

        self.resplot.set_xlim(x_min,x_max)
        self.resplot.set_ylim(y_min,y_max)
        
        self.update_plot()
    
    def update_plot(self):
        self.fig.canvas.draw()
        self.fig.tight_layout()

    def update_legend(self,label=None):

        _handles=[]

        if self.lines_data_list!=[] and self.lines_data_list is not None:
            _handles.extend(self.lines_data_list)


        if self.lines_model_list!=[] and self.lines_model_list is not None:
            _handles.extend(self.lines_model_list)



        for h in _handles[:]:
            #print('_label',h._label)
            if h._label is  None:
                _handles.remove(h)
            elif h._label.startswith('_line'):
                _handles.remove(h)
            else:
                 pass

        #for h in _handles:
        #    print('_label',h._label)

        self.sedplot.legend(handles=_handles,loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1, prop={'size':10})
        self.update_plot()




    def add_model_plot(self, model, label=None, color=None, line_style=None, flim=None,auto_label=True,fit_range=None):

        try:
            #print( "a")
            x, y = model.get_model_points(log_log=True, frame = self.frame)
        except Exception as e:
            #print("a",e)
            try:

                x, y = model.SED.get_model_points(log_log=True, frame = self.frame)

            except Exception as e:
                #print("b", e)
                #print(model, "!!! Error has no SED instance or something wrong in get_model_points()")
                #print(e)
                return

        #print('->x,y',x,y)
        #if color is None:
        #    color = self.counter

        if line_style is None:
            line_style = '-'

        if label is None and auto_label is True:
            if model.name is not None:
                label = model.name
            else:
                label = 'line %d' % self.counter
        if flim is not None:

            msk=y>np.log10(flim)
            x=x[msk]
            y=y[msk]
        else:
            pass

        if fit_range is not None:
            msk1 = x < fit_range[1]
            msk2 = x > fit_range[0]

            x = x[msk1 * msk2]
            y = y[msk1 * msk2]

        line, = self.sedplot.plot(x, y, line_style, label=label,color=color,linewidth=1.0)


        self.lines_model_list.append(line)

        self.update_legend()
        self.update_plot()

        self.counter += 1

    def add_data_plot(self,sed_data,label=None,color=None,autoscale=True,fmt='o',ms=4,mew=0.5,fit_range=None):



        try:
            x,y,dx,dy,=sed_data.get_data_points(log_log=True,frame=self.frame)
        except:
            print ("!!! ERROR failed to get data points from", sed_data)
            print
            raise RuntimeError
            
         
       
        # get x,y,dx,dy from SEDdata
        if dx is None:
            dx=np.zeros(len(sed_data.data['nu_data']))
        

        if dy is None:
            dy=np.zeros(len(sed_data.data['nu_data']))

        UL = sed_data.data['UL']
          
        # set color
        #if color is None:
        #    color=self.counter

                
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
        self.update_legend()
        self.update_plot()
        

    def add_xy_plot(self,x,y,label=None,color=None,line_style=None,autoscale=False):

        #color setting  
        #if color is None:
        #    color=self.counter
       
    
        
        if line_style is None:
            line_style='-'
           

        if label is None:
            label='line %d'%self.counter

        line, = self.sedplot.plot(x, y, line_style,label=label)

        

        self.lines_model_list.append(line)


        self.counter+=1

        self.update_legend()
        self.update_plot()



       
                
        

    def add_residual_plot(self,model,data,label=None,color=None,filter_UL=True,fit_range=None):

        if self.counter_res == 0:
            self.add_res_zeroline()
        #print('bbbbb')
        if data is not None:

            x,y=model.get_residuals(log_log=True,data=data,filter_UL=filter_UL)

            if fit_range is not None:
                msk1 = x<fit_range[1]
                msk2 =x>fit_range[0]

                x=x[msk1*msk2]
                y=y[msk1*msk2]
            #print('aaaaaaa',fit_range,x)
            line = self.resplot.errorbar(x, y, yerr=np.ones(x.size), fmt='+',color=color)
            self.lines_res_list.append(line)
            self.counter_res += 1
        else:
            pass




        




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


    def save(self,filename):
        outfile=self.out_dir+'/'+filename
        self.fig.savefig(outfile)


    def show(self):
        self.fig.show()




class  PlotPdistr (object):

    def __init__(self):
        self.fig, self.ax = plt.subplots()

    def plot_distr(self,gamma,n_gamma,y_min=None,y_max=None,x_min=None,x_max=None):

        self.ax.plot(np.log10(gamma), np.log10(n_gamma))
        self.ax.set_xlabel(r'log($\gamma$)')
        self.ax.set_ylabel(r'log(n($\gamma$))')
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlim(x_min, x_max)
        self.update_plot()

    def plot_distr3p(self,gamma,n_gamma,y_min=None,y_max=None,x_min=None,x_max=None):

        self.ax.plot(np.log10(gamma), np.log10(n_gamma *  gamma * gamma *  gamma))
        self.ax.set_xlabel(r'log($\gamma$)')
        self.ax.set_ylabel(r'log(n($\gamma$) \gamma^3)')
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlim(x_min, x_max)
        self.update_plot()

    def update_plot(self):
        self.fig.canvas.draw()


class  PlotSpecComp (object):

    def __init__(self):
        self.fig, self.ax = plt.subplots()


    def plot(self,nu,nuFnu,y_min=None,y_max=None):

        self.ax.plot(np.log10(nu), np.log10(nuFnu))
        self.ax.set_xlabel(r'log($ \nu $)  (Hz)')
        self.ax.set_ylabel(r'log($ \nu F_{\nu} $ )  (erg cm$^{-2}$  s$^{-1}$)')
        self.ax.set_ylim(y_min, y_max)
        self.update_plot()


    def update_plot(self):
        self.fig.canvas.draw()


class  PlotSeedPhotons (object):

    def __init__(self):
        self.fig, self.ax = plt.subplots()


    def plot(self,nu,nuFnu,y_min=None,y_max=None):

        self.ax.plot(np.log10(nu), np.log10(nuFnu))
        self.ax.set_xlabel(r'log($ \nu $)  (Hz)')
        self.ax.set_ylabel(r'log(n )  (photons cm$^{-3}$  Hz$^{-1}$  ster$^{-1}$)')
        self.ax.set_ylim(y_min, y_max)
        self.update_plot()


    def update_plot(self):
        self.fig.canvas.draw()