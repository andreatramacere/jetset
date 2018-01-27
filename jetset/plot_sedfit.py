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

NOPYLAB=True

try:
    
    import  pylab as plt
   
    from matplotlib import pylab as pp
except:

    NOPYLAB=True

    #print "pylab not found on this system"
    #print "install package, and/or update pythonpath"



    
from collections import namedtuple

from output import section_separator,WorkPlace

import numpy as np


__all__=['Plot']


class  Plot (object):
    
    
    def __init__(self,SEDdata=None,x_min=None,x_max=None,y_min=None,y_max=None,interactive=True,plot_workplace=None,title='Plot'):
      
        self.axis_kw=['x_min','x_max','y_min','y_max']
        self.interactive=interactive
        self.x_min=x_min
        self.x_max=x_max
        self.y_min=y_min
        self.y_max=y_max
        plot_workplace=plot_workplace
        self.line_tuple=namedtuple('line',['label','ref'])
        self.lines_data_list=[]
        self.lines_model_list=[]
        self.legend=[]
     
        
        #plt.ioff()
        if self.interactive==True:
       
            
            plt.ion()
            print 'running PyLab in interactive mode'
        

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
        
        try:
             import matplotlib.gridspec as gridspec
             GRIDSPEC=True
             gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
        
        except:
             GRIDSPEC=False
             gs=[]
             gs.append(plt.subplot(2,1,1))

             

            
            
        
        self.sedplot= self.fig.add_subplot(gs[0])
        
        
        self.set_plot_axis_labels(SEDdata)
        
       
        self.sedplot.set_autoscale_on(False)
        
        self.counter=0
        

        self.sedplot.grid(True)   
       
        #Build residuals plot
        self.x_min_res=self.x_min
        self.x_max_res=self.x_max
        
        self.y_min_res=-1
        self.y_max_res=1
        
            
        if SEDdata is not None :
            if GRIDSPEC==True:
                self.resplot= self.fig.add_subplot(gs[1])
            else:
                 gs.append(plt.subplot(2,1,2))
                 self.resplot= self.fig.add_subplot(gs[1])

            
            self.lx_res='log($ \\nu $)  (Hz)'    
            self.ly_res='res'
        
            
            self.resplot.set_ylabel(self.ly_res)
            self.resplot.set_xlabel(self.lx_res)
            
            self.resplot.set_xlim(self.x_min_res,self.x_max_res)
            self.resplot.set_ylim(self.y_min_res,self.y_max_res)
            
            self.add_res_zeroline()
            
            
        else:
            self.resplot=None
            
        self.counter_res=0
        #self.fig.show()
    
    
    def clean_data_lines(self):
        
        for i in xrange(len(self.lines_data_list)):
            self.del_data_line(0)
    
    def clean_model_lines(self):
        for i in xrange(len(self.lines_model_list)):
            self.del_model_line(0)
            
            
    def list_lines(self):
        
        if self.lines_data_list==[] and self.lines_model_list==[]:
            print  "no lines to show "
        
        else:
            ID=0
            print "data"
            for plot_line in self.lines_data_list:
                
                print ID, plot_line.label, plot_line.ref
                
                
                ID+=1
            
            ID=0    
            print "models"
            for plot_line in self.lines_model_list:
                
                print ID, plot_line.label, plot_line.ref
                
                
                ID+=1
    
    def del_data_line(self,line_ID):
        
        if self.lines_data_list==[]:
            print  "no lines to delete "
        
        else:
            
            print "removing line: ",self.lines_data_list[line_ID]
            
        line=self.lines_data_list[line_ID][1]
        
        for item in line:
            #This removes lines
            if np.shape(item)==():
                item.remove()
            else:
                #This removes containers for data with errorbars
                for item1 in item:
                    item1.remove()
            
            
            #self.sedplot.lines.remove(self.lines_list[line_ID].ref[0])
                    
        self.legend.remove(self.lines_data_list[line_ID].label)
 
        del self.lines_data_list[line_ID]
 
        #print self.legend
        
        self.update_plot()
    
    
    
        
    def update_plot(self):
     

       self.fig.canvas.draw()
       #self.fig.canvas.draw()
       #self.fig.show()     
       #plt.ioff()
    
    
    
    def del_model_line(self,line_ID):
        
        if self.lines_model_list==[]:
            print  "no lines to delete "
        
        else:
            
            print "removing line: ",self.lines_model_list[line_ID]
            
        line=self.lines_model_list[line_ID][1]
        
        for item in line:
            #This removes lines
            if np.shape(item)==():
                item.remove()
            else:
                #This removes containers for data with errorbars
                for item1 in item:
                    item1.remove()
            
            
            #self.sedplot.lines.remove(self.lines_list[line_ID].ref[0])
                    
        self.legend.remove(self.lines_model_list[line_ID].label)
 
        del self.lines_model_list[line_ID]
 
        #print self.legend
        
        self.update_legend()
        self.update_plot()
        
        
    
    
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
        
        self.sedplot.set_xlim(self.x_min,self.x_max)
        self.sedplot.set_ylim(self.y_min,self.y_max)
        
    
    def add_res_zeroline(self):
        
        y0=np.zeros(2)
        x0=[self.x_min_res,self.x_max_res]
       
       
        self.resplot.plot(x0,y0,'--',color='black')
       
        self.update_plot()
        
        
        
     
     
     
     
    def rescale(self,**kw):

        """

        Rescales the data plot acording to x_min/max, y_min/max values
        if no argument are passed, the plot is rescaled according
        to current x_min/max, y_min/max values
        if a value is passed as argument, the corresponding axis 
        min/max is updated, and the plot is rescaled accordingly
        
        Args: x_min,x_max,y_min,y_max

        """
        
        
        #  check if user did provide
        #  keywords==self.axis_kw
        #  and updates accordingly
        #  the following min/max values 
        #  are left unchanged
 
        keys=kw.keys()
        
        for k in keys:
            
            if k in self.axis_kw:
                setattr(self,k,kw[k])
            else:
                print "the keyword %s is not in allowed=%s"%(k,self.axis_kw)
                print "please correct"  
                print 

                return 
            

            
        self.sedplot.set_xlim(self.x_min,self.x_max)
        self.sedplot.set_ylim(self.y_min,self.y_max)
        
        self.update_plot()
        
    
    def autoscale(self):
        self.sedplot.set_autoscale_on(True)
        
        self.sedplot.autoscale_view(tight=True)
        
        self.x_min,self.x_max=self.sedplot.get_xlim()
        
        self.y_min,self.y_max=self.sedplot.get_ylim()
        
        self.sedplot.set_xticks(np.arange(int(self.x_min)-2,int(self.x_max)+2,1.0))
        
        self.sedplot.set_xlim(self.x_min-1,self.x_max+1)
        
        self.sedplot.set_ylim(self.y_min-1,self.y_max+1)
        
        
        
        if self.resplot is not None  :
            
            #self.resplot.autoscale_view(tight=True)
            
            self.x_min_res=self.x_min-1
            self.x_max_res=self.x_max+1

            self.y_min_res,self.y_max_res=self.resplot.get_ylim()
            
            self.resplot.set_xticks(np.arange(int(self.x_min_res)-2,int(self.x_max_res)+2,1.0))
            
            self.resplot.set_xlim(self.x_min_res,self.x_max_res)
            self.resplot.set_ylim(self.y_min_res,self.y_max_res)
            
            
            
        self.update_plot()
        
        self.sedplot.set_autoscale_on(False)
   
   
    
    def rescale_res(self,**kw):
    
        """

        Rescales the data plot acording to x_min/max, y_min/max values
        if no argument are passed, the plot is rescaled according
        to current x_min/max, y_min/max values
        if a value is passed as argument, the corresponding axis 
        min/max is updated, and the plot is rescaled accordingly
        
        Args: x_min,x_max,y_min,y_max

        """
        
        
        #  check if user did provide
        #  keywords==self.axis_kw
        #  and updates accordingly
        #  the following min/max values 
        #  are left unchanged
 
        keys=kw.keys()
        
        for k in keys:
            
            if k in self.axis_kw:
                setattr(self,k,kw[k])
            else:
                print "the keyword %s is not in allowed=%s"%(k,self.axis_kw)
                print "please correct"  
                print 

                return 
            

            
        self.resplot.set_xlim(self.x_min_res,self.x_max_res)
        self.resplot.set_ylim(self.y_min_res,self.y_max_res)
        
        self.update_plot()
    
            
    def update_legend(self,label=None):
        
        if label is not None  :
            
            self.legend.append(label)


        self.sedplot.legend(self.legend,loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3, prop={'size':12})
    
        #self.fig.show()

        #self.update_plot()
    
        
        
    def add_data_plot(self,SEDdata,label=None,color=None,autoscale=False):
        try:
            x,y,dx,dy,=SEDdata.get_data_points(log_log=True)
        except:
            print "!!! ERROR failed to get data points from", SEDdata
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
        
        
        self.sedplot.set_autoscale_on(False)
        
        line=self.sedplot.errorbar(x, y,xerr=dx, yerr=dy, fmt='o',lolims=SEDdata.data['UL'])
        

       
        
                
        if label is None:
            if SEDdata.obj_name is not None  :
                label=SEDdata.obj_name
            else:
                label='line %d'%self.counter
        
        self.update_legend(label)
        
        newline=self.line_tuple(label,line)
        
        self.lines_data_list.append(newline)
        
        if autoscale==True:
            self.autoscale()
        
        self.counter+=1
        
        self.update_plot()
        

    def add_xy_plot(self,x,y,label=None,color=None,line_style=None,autoscale=False):
       
       
        #color setting  
        if color is None:
            color=self.counter
       
    
        
        if line_style is None:
            line_style='-'
           
        line=self.sedplot.plot(x,y,line_style)
       
        line[0].set_ydata(y)
        line[0].set_xdata(x)
         
        if label is None:
            label='line %d'%self.counter

        self.update_legend(label)
        
        newline=self.line_tuple(label,line)
        
        self.lines_model_list.append(newline)

        if autoscale==True:
            self.autoscale()
        
        self.counter+=1
    
        self.update_plot()       
        
        
        
        
          
    def add_model_plot(self,model,label=None,color=None,line_style=None,autoscale=False,update=True):
    
        try:
            #print "a"
            x,y=model.get_model_points(log_log=True)
        except:
            try:
                #print "b"
                x,y=model.SED.get_model_points(log_log=True)
            except:
                print model, "!!! Error has no SED instance or something wrong in get_model_points()"
                return
        
        #print "x,y",model,x,y
        #color setting  
        if color is None:
            color=self.counter
       
       
       
       
        if line_style is None:
            line_style='-'
        
        self.sedplot.set_autoscale_on(False)
         
        line=self.sedplot.plot(x,y,line_style)
        #print line
        line[0].set_ydata(y)
        line[0].set_xdata(x)
        
            
        if label is None:
            if model.name is not None  :
                label=model.name
            else:
                label='line %d'%self.counter
        
        
        
        newline=self.line_tuple(label,line)
        
        self.lines_model_list.append(newline)
        
        self.update_legend(label)
          
        if autoscale==True:
            self.autoscale()
        
        self.counter+=1
        
        
        if update==True:
            self.update_plot()
            
       
                
        

    def add_residual_plot(self,model,label=None,color=None,autoscale=False):
        
        try:
            x,y=model.get_residuals(log_log=True)
        except:
            try:
                x,y=model.SED.get_residuals(log_log=True)
            except:
                print model, "has no residuals"
                return

      
        if self.counter_res==0:
            self.add_res_zeroline()
        
        self.resplot.set_autoscaley_on(True)

        self.resplot.plot(x,y,'+')
           
        
        
        if autoscale==True:
            self.resplot.autoscale()
    
        self.counter_res+=1

        self.update_plot()
          
        
    
    
    
    def add_text(self,lines):
        self.PLT.focus(0,0)

        t=''
        for line in lines:
            t+='%s \\n'%line.strip()
        self.PLT.text(t,font=10,charsize=0.6,x=self.x_min-1.5,y=self.y_min-2.85)
        self.PLT.redraw()


    def save(self,filename):
        outfile=self.out_dir+'/'+filename
        self.fig.savefig(outfile)
