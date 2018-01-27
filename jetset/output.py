#!/usr/bin/env python
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


__all__=['clean_dir','makedir','section_separator', 'WorkPlace']




def clean_dir(dir_name):
    print "cleaning dir",dir_name
    shutil.rmtree(dir_name)
    os.mkdir(dir_name)

def makedir(out_dir,clean_work_dir=True):
    """
    creates a directory
    """
    if os.path.isdir(out_dir):
        print "directory %s already existing"%(out_dir)
        if clean_work_dir==True:
            print 'removing existing dir'
            shutil.rmtree(out_dir)
            os.mkdir(out_dir)
            print 'the directory %s has been created' % (out_dir)
    else:
        if os.path.isfile(out_dir):
            print "a file with the same name of dir=%s, exists"%out_dir
            print "select a differn name"
        else:
            os.mkdir(out_dir)
            print 'the directory %s has been created'%(out_dir)


class WorkPlace(object):
    """
    Class to set the working place
    using static members (class members)
    
    Variables
    
    :ivar out_dir: directory name (default=./) 
    :ivar flag: flag name (default=sed-fit-test)
    
    
    """
    def __init__(self):
        self.out_dir='./'
        self.flag='sed-fit-test'
    
    
    
    
            
            

section_separator=    "=============================================================================================\n"
