"""
===================================================================
Moudule: plot_sedfit
===================================================================

This module contains all the classes necessary to build 


..



Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.plot_sedfit
    


  
.. autosummary::
  
   
   
    
Module API
-------------------------------------------------------------------

"""


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

import numpy as np


__all__=['Plot']


class  Plot (object):
    
    
    def __init__(self,SEDdata=None,
                 model=None,
                 x_min=6.0,
                 x_max=3.0,
                 y_min=-20.0,
                 y_max=-9.0,
                 interactive=False,
                 plot_workplace=None,
                 title='Plot'):
                 #autoscale=True):
      
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
      
            
        self.fig=plt.figure(figsize=(12,9))
        


        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
            
        
        self.sedplot= self.fig.add_subplot(gs[0])
        
        
        self.set_plot_axis_labels(SEDdata)
        
        #if autoscale==True:
        self.sedplot.set_autoscalex_on(True)
        self.sedplot.set_autoscaley_on(True)
        self.sedplot.set_autoscale_on(True)
        self.counter=0
        

        self.sedplot.grid(True)   
       
        #Build residuals plot

        
            
        if SEDdata is not None :
            self.add_data_plot(SEDdata)

        if model is not  None:
            self.add_model_plot(model)

        self._add_res_plot()
        self.resplot.set_autoscaley_on(True)
            
        #else:
        #    self.resplot=None
            
        self.counter_res=0
        #self.fig.show()
        #if autoscale==True:
        #    self.autoscale()


    def _add_res_plot(self):



        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
        self.resplot = self.fig.add_subplot(gs[1], sharex=self.sedplot)

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

    def set_plot_axis_labels(self,SEDdata=None):
        
        if SEDdata is not None :
            self.lx='log($ \\nu $)  (Hz)'
                
            if SEDdata.restframe=='src':
                self.ly='log($ \\nu L_{\\nu} $ )  (erg  s$^{-1}$)' 
            
            elif SEDdata.restframe=='obs':
                self.ly='log($ \\nu F_{\\nu} $ )  (erg cm$^{-2}$  s$^{-1}$)' 
        
        else:
            self.lx='log($ \\nu $)  (Hz)'
            self.ly='log($ \\nu F_{\\nu} $ )  (erg cm$^{-2}$  s$^{-1}$)'
            
        self.sedplot.set_ylabel(self.ly)
        self.sedplot.set_xlabel(self.lx)
        
        #self.sedplot.set_xlim(self.x_min,self.x_max)
        #self.sedplot.set_ylim(self.y_min,self.y_max)
        
    
    def add_res_zeroline(self):
        x_min_res, x_max_res = self.resplot.get_xlim()


        y0=np.zeros(2)
        x0=[x_min_res,x_max_res]
       
        print(x0,y0)
        self.resplot.plot(x0,y0,'--',color='black')
        self.update_plot()

        
        
     
     
     
     
    def rescale(self,x_min=None,x_max=None,y_min=None,y_max=None):


            
        self.sedplot.set_xlim(x_min,x_max)
        self.sedplot.set_ylim(y_min,y_max)


    
    def autoscale(self):

        
        #self.sedplot.autoscale_view(tight=True)
        for l in self.sedplot.lines:

            x_min,x_max=self.sedplot.get_xlim()

            y_min,y_max=self.sedplot.get_ylim()
        
        self.sedplot.set_xticks(np.arange(int(x_min)-2,int(x_max)+2,1.0))
        
        self.sedplot.set_xlim(x_min-1,x_max+1)
        
        self.sedplot.set_ylim(y_min-1,y_max+1)
        
        
        
        if self.resplot is not None  :
            
            #self.resplot.autoscale_view(tight=True)
            
            #self.x_min_res=self.x_min-1
            #self.x_max_res=self.x_max+1

            y_min_res,y_max_res=self.resplot.get_ylim()
            x_min_res,x_max_res=self.resplot.get_xlim()

            
            self.resplot.set_xticks(np.arange(int(x_min_res)-2,int(x_max_res)+2,1.0))
            
            self.resplot.set_xlim(x_min_res,x_max_res)
            self.resplot.set_ylim(y_min_res,y_max_res)
            
            
            
        #self.update_plot()
        
        #self.sedplot.set_autoscale_on(False)
   
   
    
    def rescale_res(self,x_min=None,x_max=None,y_min=None,y_max=None):

        self.resplot.set_xlim(x_min,x_max)
        self.resplot.set_ylim(y_min,y_max)
        
        self.update_plot()
    
    def update_plot(self):
        self.fig.canvas.draw()

    def update_legend(self,label=None):
        

        _handles=[]
        _labels=[]
        if self.lines_data_list!=[] and self.lines_data_list is not None:
            _handles.extend(self.lines_data_list)

        if self.lines_model_list!=[] and self.lines_model_list is not None:
            _handles.extend(self.lines_model_list)


        if _handles==[]:
            _handles=None

        if _handles is None:
            _labels=[]
        else:
            _labels=None

        #print('_handles', _handles, _labels)
        self.sedplot.legend(handles=_handles,labels=_labels,loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3, prop={'size':12})
        self.update_plot()




    def add_model_plot(self, model, label=None, color=None, line_style=None, autoscale=False, update=True):

        try:
            # print "a"
            x, y = model.get_model_points(log_log=True)
        except:
            try:
                # print "b"
                x, y = model.SED.get_model_points(log_log=True)
            except:
                print(model, "!!! Error has no SED instance or something wrong in get_model_points()")
                return

        if color is None:
            color = self.counter

        if line_style is None:
            line_style = '-'

        if label is None:
            if model.name is not None:
                label = model.name
            else:
                label = 'line %d' % self.counter

        line, = self.sedplot.plot(x, y, line_style, label=label)


        self.lines_model_list.append(line)

        self.update_legend()
        self.update_plot()

        self.counter += 1

    def add_data_plot(self,SEDdata,label=None,color=None,autoscale=True,fmt='o',ms=None,mew=None):



        try:
            x,y,dx,dy,=SEDdata.get_data_points(log_log=True)
        except:
            print ("!!! ERROR failed to get data points from", SEDdata)
            print
            raise RuntimeError
            
         
       
        # get x,y,dx,dy from SEDdata
        if dx is None:
            dx=np.zeros(len(SEDdata.data['nu_data']))
        

        if dy is None:
            dy=np.zeros(len(SEDdata.data['nu_data']))
        
        
          
        # set color
        if color is None:
            color=self.counter

                
        if label is None:
            if SEDdata.obj_name is not None  :
                label=SEDdata.obj_name
            else:
                label='line %d'%self.counter

        line = self.sedplot.errorbar(x, y, xerr=dx, yerr=dy, fmt=fmt
                                     , uplims=SEDdata.data['UL'],label=label,ms=ms,mew=mew)


        


        self.lines_data_list.append(line)



        self.counter+=1
        self.update_legend()
        self.update_plot()
        

    def add_xy_plot(self,x,y,label=None,color=None,line_style=None,autoscale=False):

        #color setting  
        if color is None:
            color=self.counter
       
    
        
        if line_style is None:
            line_style='-'
           

        if label is None:
            label='line %d'%self.counter

        line, = self.sedplot.plot(x, y, line_style,label=label)

        

        self.lines_model_list.append(line)


        self.counter+=1

        self.update_legend()
        self.update_plot()



       
                
        

    def add_residual_plot(self,model,data,label=None,color=None,filter_UL=True):



        if data is not None:
            x,y=model.get_residuals(log_log=True,data=data,filter_UL=filter_UL)
        else:
            pass


        if self.counter_res==0:
            self.add_res_zeroline()

        

        line = self.resplot.errorbar(x, y, yerr=np.ones(x.size), fmt='+')

        self.lines_res_list.append(line)
        self.update_plot()
        self.counter_res+=1





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
