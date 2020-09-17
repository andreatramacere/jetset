#from __future__ import absolute_import, division, print_function

#from builtins import (bytes, str, open, super, range,
#                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

import os
import sys
import glob

from .cosmo_tools import  Cosmo
from string import *
import math as m
from astropy.io import fits as pyfits
import numpy as n

__all__=['check_is_number','do_convert']

def check_is_number(t):
    try: 
        float(t[0])
        return True
    except ValueError:
        return False



#------ Main --------------
def main(argv=None,cmd_par=None):

    if cmd_par  is None:
        argc=len(sys.argv)
        arg_list=sys.argv
        in_seds=arg_list[1:]
    else:
        argc=len(cmd_par)
        arg_list=(cmd_par)
        if '*' in arg_list[1]:
            in_seds=glob.glob(arg_list[1])

        print (in_seds)
 
    print (arg_list)
    if argc==1:
        print("./concvert_qdp.py infiles.qdq")
        print("exmaple: ./concvert_qdp.py SED_Ciprini_v1/sed-*.qdp")
        print("input data are in obs rest frame, out data in src rest frame (z-corrected)")
        sys.exit(0)

   
    print (arg_list)

    comos_eval=Cosmo(units='cm')
    for sed in in_seds:
        do_convert(sed,comos_eval)
        
def do_convert(infile_name,cosmo_eval):
    infile=open(infile_name,'r')
    
    outfile_name=infile_name.split('.')
    outfile_name[-1]='dat'
    outfile_name=".".join(outfile_name)
    outfile=open(outfile_name,'w')
    obj_name=infile_name.split('/')[-1]
    obj_name=obj_name.split('.')[0]
    obj_name=obj_name.split('-')[1]
    
    inlines=infile.readlines()
    
    infile.close()
    
    print('# set_number: -1 hyst, 0 sym, 1 green, 2 blue ',file=outfile)
    print('# file_name ',obj_name,file=outfile)
    set_code={}
    UL_sets=[]

    for line in inlines:
        line=line.strip()
        if len(line)>0:
            tkn=line.split()
            #print tkn
            if tkn[0]=="!" and tkn[1]=='Redshift':
                z=float(tkn[3])
                DL=cosmo_eval.DL(z)
                print('# z ',z,file=outfile)
                print('# DL(cm) ',DL,file=outfile)
                print('# restframe  src',file=outfile)
                print('# dataScale  log-log',file=outfile)
            
            for t in tkn:
                #print t
                if ('J' in t and '+' in t) or ('J' in t  and '-' in t):
                    name=t
                    #print 'J' and '+' in t
                    #print 'J' and '-' in t 
                    #print t
                    print('# obj_name    ',name,file=outfile)
    
            if tkn[0]=='col' and tkn[2]=='on':
                #print tkn,tkn[3:]
                for  t in tkn[3:]:
                    #print t,tkn[1],t=='2'
                    #sym
                    if tkn[1]== '2':
                        set_code[t]='0'
                        
                    #
                    elif tkn[1]=='15':
                        set_code[t]='-1'
                    #green
                    elif tkn[1]=='3':
                        set_code[t]='1'
                    #blue
                    elif tkn[1]=='4':
                        set_code[t]='2'

            #UL
            if tkn[0]=='ma' and tkn[1]=='31':
                UL_sets=tkn[3:]
                print ("UL sets",UL_sets)


    print>>outfile,'#nu nuLnu dnu dnuLnu set_number '     
    set=1
    #print set_code
    if z!=0.0:
        for line in inlines:
            print (tkn)
            tkn=line.split()           
            if  len(tkn)==2 and tkn[0]=='NO' and tkn[1]=='NO':
                #print tkn
                set+=1
                #print>>outfile,"###############################"
            if check_is_number(tkn[0]):
                if float(tkn[0])>5:
                    #print tkn
                    nu_obs=pow(10,float(tkn[0]))
                    dnu_obs=nu_obs*float(tkn[1])
                    nuLnu_obs=pow(10,float(tkn[2]))
                
                    
                    nu_rest=nu_obs*(1+z)
                    log_dnuLnu_rest=float(tkn[3])
                    #!!! ASSUMES UL in FERMI if dy==0
                    #!!! for other data sets UL from set with ma=31
                    if nu_obs>1.0e22 and log_dnuLnu_rest==0.0:
                        #print  dnu_rest
                        log_dnuLnu_rest=-1
                    
                    for UL_set in UL_sets:
                        if int(UL_set)==set:
                            #print 'UL',UL_set,set
                            log_dnuLnu_rest=-1

                    nuLnu_rest=nuLnu_obs*4*m.pi*DL*DL
                    #errore relativo 
                    #dnuLnu_rest=nuLnu_rest*(float(tkn[3]))
                        
                    print(m.log10(nu_rest),m.log10(nuLnu_rest),0,log_dnuLnu_rest,set_code['%s'%set],file=outfile)
                    #print outfile,nu_rest,dnu_rest,nuLnu_rest,dnuLnu_rest
            
                
        outfile.close()

if __name__ == "__main__":
    main()
