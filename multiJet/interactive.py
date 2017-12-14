#!/usr/bin/env python 

'''
Created on 2013 1 16

@author: orion
'''

from test_data_helper import test_SEDs
from data_loader import ObsData
from minimizer import fit_SED
from model_manager import FitModel
from obs_constrain import ObsConstrain
from plot_sedfit import Plot
from sed_shaper import SEDShape
from output import workplace
from cosmo_tools import Cosmo
from jet_model import Jet
from template_model import Template

 




welcome="""
***************************************************************
* Welcome to the interactive version of BlazarSEDFit          *
*                                                             *
* author andrea tramacere 2013                                *
*                                                             *
* to access the commands:                                     *
*                                                             *
*    commands(SEDFit)                                         *
***************************************************************
"""


print welcome
    
#import_dic={}
#import_dic['data_loader']=['ObsData']
#import_dic['minimizer']=['fit_par','fit_par_array','fit_SED']
#import_dic['model']=['SEDmodel']
#import_dic['obs_constrain']=['SED_obspar']
#import_dic['plot_sedfit']=['Plot']
#import_dic['sed_shaper']=['SEDShape']
#import_dic['output']=['set_workplace']
#
#for key in import_dic.keys():
#    print key,import_dic[key]  
#    __import__('BlazarSEDFit.'+key, fromlist=import_dic[key])

#class setup_interactive_session(object):
#    
#    
#    
#    import_dic={}
#    
#    import_dic['data_loader']=['ObsData']
#    
#    
#    import_dic['minimizer']=['fit_par','fit_par_array','fit_SED']
#    
#    import_dic['model']=['SEDmodel']
#
#    import_dic['obs_constrain']=['SED_obspar']
#
#    import_dic['plot_sedfit']=['Plot']
#
#    import_dic['sed_shaper']=['SEDShape']
#    
#    import_dic['output']=['set_workplace']
#                
#    
#
#
#    welcome="""
#***************************************************************
#* Welcome to the interactive version of BlazarSEDFit          *
#*                                                             *
#* author andrea tramacere 2013                                *
#*                                                             *
#* to access the commands:                                     *
#*                                                             *
#*    SEDFit.help()                                              *
#***************************************************************
#"""
#    for key in import_dic.keys():
#        print key,import_dic[key]  
#        mod=key
#        print mod
#        __import__(mod, fromlist=import_dic[key],level=1)
#        
#        
#        
#        
#        def help(self):
#            print "help"
#            for key in self.import_dic.keys():
#                print self.import_dic[key]
#            
#        
#        
        