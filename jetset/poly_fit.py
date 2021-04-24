
__author__ = "Andrea Tramacere"

from scipy.optimize import leastsq
import scipy as s
from numpy import polyfit,polyval,polyder
import  numpy as np

__all__=['check_maxima','cubic_peaks','do_cubic_fit','do_cubic_fit_peaks','do_linear_fit',
         'do_log_Parab_FIT','filter_interval','find_max_cubic','lin','parab','residuals_cubic_peaks',
         'residuals_linfit','residuals_parab']


def filter_interval(x,x_range):
    msk1=x>=x_range[0]
    msk2=x<=x_range[1]
    return msk1*msk2




#-------- LINEAR FIT----------------------------
def do_linear_fit(x,y,dy=None,x_range=None):
    
    if x_range is not None  :
        msk = filter_interval(x,x_range)
        x_fit=x[msk]
        y_fit=y[msk]
        dy_fit=dy[msk] 
    else:
        x_fit=x
        y_fit=y 
        dy_fit=dy
    
    if dy is None:
        dy_fit=s.ones(len(x_fit))
   
        


    if len(x_fit)<2:
        raise RuntimeError("fit failed number of points less then 2")
     
    #print nu_fit,nuFnu_fit,err_nuFnu_fit
    pinit = [1.0,10.0]
    #print "input", log_nu,log_nuLnu    
    p,cov,infodict,mesg,ier = leastsq(residuals_linfit, pinit,args=(x_fit,y_fit,dy_fit), full_output=1)
    
    chisq=(infodict['fvec']**2).sum()
    dof=len(x_fit)-len(p)
    if dof>0:
        rmse=s.sqrt(chisq/dof)
    else:
        rmse=0
        
    if cov is not None  :
        err=s.sqrt(cov.diagonal())*rmse
    else:
        err=s.zeros(len(p))
    #print p
    return p,err

def lin(p,x):

    return p[0]+p[1]*(x)

def residuals_linfit(p,x,y,dy):
    #print len(dy),len(x),len(y),len(p)
    err = (y-lin(p,x))/dy
    return err
#-----------------------------------------------


#----CUBIC-PEAKS------------------------------------
def cubic_peaks(p,x):
    #p0=yp
    #p1=b
    #p2=xp
    #
    return p[0] + p[1]*(x-p[2])*(x-p[2]) + p[3]*(x-p[2])*(x-p[2])*(x-p[2])


def residuals_cubic_peaks(p,x,y,dy,template):
    if template is not None  :
        model1=s.power(10,cubic_peaks(p[:4],x))
        model2=s.power(10,template.func(x,nuFnu_scale=p[4]))
        model=s.log10(model1+model2)
    else:
        model=cubic_peaks(p[:4],x)
    
    
    
    err = (y-model)/dy
    return err

def do_cubic_fit_peaks(x,y,xp,yp,dy=None,x_range=None,template=None):
    
    if x_range is not None  :
        msk = filter_interval(x,x_range)
        x_fit=x[msk]
        y_fit=y[msk]
        dy_fit=dy[msk] 
    else:
        x_fit=x
        y_fit=y 
        dy_fit=dy
    
    if dy is None:
        dy_fit=s.ones(len(x_fit))
    
    
    if template is not None  :
        pinit = [yp,0,xp,0,yp]
    else:
        pinit = [yp,0,xp,0]
    
    
    p,cov,infodict,mesg,ier = leastsq(residuals_cubic_peaks, pinit,args=(x_fit,y_fit,dy_fit,template), full_output=1)
    
    #print"cubic fit", out[0]
    chisq=(infodict['fvec']**2).sum()
    # dof is degrees of freedom
    dof=len(x_fit)-len(p)
    if dof>0:
        rmse=s.sqrt(chisq/dof)
    else:
        rmse=0
        
    #print "--> chi_red", rmse
    if cov is not None  :
        err=s.sqrt(cov.diagonal())*rmse
    else:
        err=s.zeros(len(p))
    #print p[3]
    #print "--> fit log: ",mesg,ier
    #print infodict['fvec']
    par_parab=[p[2],p[0],p[1]]
    err_parab=[err[2],err[0],err[1]]
    if len(p)==5:
        par_parab.append(p[4])
        err_parab.append(err[4])
        
    return par_parab,err_parab,dof
#------------------------------------------------------


#-----CUBIC   -----------------------------------------
def find_max_cubic(p,x,y,dy=None,x_range=None):
    der=polyder(p,1)
    delta=der[1]*der[1] -(4*der[0]*der[2])
    if delta>=0:
        #print p,der
        root1=(-der[1]+s.sqrt(delta))/(2*der[0])
        root2=(-der[1]-s.sqrt(delta))/(2*der[0])
        #print root1,root2
        if der[0]>0:
            xp=min(root1,root2)
        else:
            xp=max(root1,root2)
        
        yp=polyval(p,xp)
        #p0=yp
        #p1=xp
        #p2=b
        #print "xp,yp=",xp,yp
        p_peaks,err,dof=do_cubic_fit_peaks(x,y,xp,yp,dy=dy,x_range=x_range)
        return  p_peaks,err,dof
    else:
        print ("!!! no maxima found for cubic fit")
        xp=None
        yp=None
        b=None
        return [yp,xp,b],[0,0,0],'No'


def check_maxima(xm,ym,b,x,y):
    c = np.array([x, y]).T
    c = c[np.argsort(c[:, 1])]
    max=c[-1]
    #print "check max", c
    #print max[0],max[1]
    #min=c[0]


    if b<0:
       return max[0],max[1]
    if b>=0:
        if ym<max[1]:
            return max[0],max[1]
        else:
            return xm,ym


def do_cubic_fit(x,y,dy=None,x_range=None):
    
    if x_range is not None  :
        msk = filter_interval(x,x_range)
        x_fit=x[msk]
        y_fit=y[msk]
        dy_fit=dy[msk] 
    else:
        x_fit=x
        y_fit=y 
        dy_fit=dy
    
    if dy is None:
        dy_fit=s.ones(len(x_fit))
    
    
    
    
    p=polyfit(x_fit,y_fit,3)
        
    val,err,dof=find_max_cubic(p,x_fit,y_fit,dy=dy_fit)
    return p,val,err
#-----------------------------------------------

#---LOG PAR-----------------------------------------

def do_log_Parab_FIT(x,y,xp,yp,beta,dy=None,x_range=None,template=None):

    if x_range is not None  :
        msk = filter_interval(x,x_range)
        x_fit=x[msk]
        y_fit=y[msk]
        dy_fit=dy[msk] 
    else:
        x_fit=x
        y_fit=y 
        dy_fit=dy
    
    if dy is None:
        dy_fit=s.ones(len(x_fit))
    
    
  
    
    if template is not None  :
        pinit = [xp,yp,beta,yp]
    else:
        pinit = [xp,yp,beta]
 
    p,cov,infodict,mesg,ier = leastsq(residuals_parab, pinit,args=(x_fit, y_fit,dy_fit,template), full_output=1)
     
    chisq=(infodict['fvec']**2).sum()

    dof=len(x_fit)-len(p)
    
    if dof>0:
        rmse=s.sqrt(chisq/dof)
    else:
        rmse=0
        
    
    if cov is not None  :
        err=s.sqrt(cov.diagonal())*rmse
    else:
        err=s.zeros(len(p))
    
    
    #print p,err
    return p,err

def parab(p,x):
    
    return p[1]+p[2]*(x-p[0])*(x-p[0])


def residuals_parab(p,x,y,dy,template):
    #print p
    if template is not None  :
        model1=s.power(10,parab(p[:3],x))
        model2=s.power(10,template.interp_template(template.func,x,p[3]))
        model=s.log10(model1+model2)
    else:
        model=parab(p,x)

    err = (y-model)/dy
    return err

#---------------------------------------------------