from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

"""
===================================================================
Moudule: output
===================================================================

This module contains all the classes necessary to the output



Classes and Inheritance Structure
-------------------------------------------------------------------

.. inheritance-diagram:: BlazarSEDFit.output
   


Classes relations
----------------------------------------------

.. figure::  classes_output.png    
   :align:   center     


  
.. autosummary::
   
   
    
Module API
-------------------------------------------------------------------
"""

import os
import shutil



__author__ = "Andrea Tramacere"

__all__=['clean_dir','makedir','section_separator', 'WorkPlace']




def clean_dir(dir_name):
    print ("cleaning dir",dir_name)
    shutil.rmtree(dir_name)
    os.mkdir(dir_name)

def makedir(out_dir,clean_work_dir=True):
    """
    creates a directory
    """
    if os.path.isdir(out_dir):
        Warning ("directory %s already existing"%(out_dir))
        if clean_work_dir==True:
            Warning ('removing existing dir')
            shutil.rmtree(out_dir)
            os.mkdir(out_dir)
            Warning ('the directory %s has been created' % (out_dir))
    else:
        if os.path.isfile(out_dir):
            Warning ("a file with the same name of dir=%s, exists"%out_dir)
            Warning ("select a differn name")
        else:
            os.mkdir(out_dir)
            Warning ('the directory %s has been created'%(out_dir))


class WorkPlace(object):
    """
    Class to set the working place
    using static members (class members)
    
    Variables
    
    :ivar out_dir: directory name (default=./) 
    :ivar flag: flag name (default=sed-fit-test)
    
    
    """
    def __init__(self,out_dir='./',flag='sed-fit-test',clean=False):
        self.out_dir=out_dir
        self.flag=flag


        makedir(out_dir,clean_work_dir=clean)

    
    
            
            

section_separator=    "===================================================================================================================\n"
