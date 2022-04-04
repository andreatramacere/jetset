
__author__ = "Andrea Tramacere"



import numpy as np

from astropy.table import Table

from .template_model import SpectralTemplateLogLog
from .output import section_separator
from .model_manager import FitModel
from .loglog_poly_model import LogParabolaEp,LogCubic,find_max_cubic,LogLinear
from .analytical_model import Disk
from .minimizer import fit_SED
from .plot_sedfit import  PlotSED

__all__=['filter_interval','find_E0','IC_fit_range','index','index_array','index_typecasting','peak_values',
         'SEDShape','spectral_index_range','sync_fit_range']

def filter_interval(x,x_range):
    msk1=x>=x_range[0]
    msk2=x<=x_range[1]
    return msk1*msk2



#----------------------------------------------------
class index_typecasting(object):
    """
    Class to handle different types of spectral indices
    """
    def __init__(self,data_type,val,error=False):
        self.spectral=None
        self.photon=None
        self.sed=None
        
         
        if error==False:
            if data_type=='spectral':
                self.spectral=val
                self.photon=val-1
                self.sed=val+1
            
            if data_type=='photon':
                self.spectral=val+1
                self.photon=val
                self.sed=val+2
    
            if data_type=='sed':
                self.spectral=val-1
                self.photon=val-2
                self.sed=val
        else:
            self.spectral=val
            self.photon=val
            self.sed=val
    

class index(object):
    """
    Class for the spectral indices
    """
    
    def __init__(self,name=None,data_type=None,val=None,err=None,idx_range=[]):
        index_names=['radio','radio_mm','mm_IR','IR_Opt','Opt_UV','BBB','X','UV_X','Fermi','TeV']

        data_type_allowed=['spectral','photon','sed']
        
        if name in index_names:
            self.name=name
        else:
            raise RuntimeError("index name=%s not allowed, allowed names=%"%(name,index_names))

        if data_type in data_type_allowed:
            self.data_type=data_type
        else:
            raise RuntimeError("idx_type =%s not allowed, possible values are=%s"%(data_type,data_type_allowed))

        if val is not None  :
            self.val=index_typecasting(self.data_type,val=val)
        else:
            self.val=None
        
        if err is not None  :
            self.err=index_typecasting(self.data_type,val=err,error=True)
        else:
            self.err=None
        
        if idx_range!=[]:
            self.idx_range=tuple(idx_range)
        else:
            self.idx_range=spectral_index_range(name)
            
    
    def assign_val(self,val,err):
        self.val=index_typecasting(self.data_type,val=val)
        self.err=index_typecasting(self.data_type,val=err,error=True)
       
    def show_val(self):
        if self.val is None:
            val='No'
        else:
            val=str("%e"%self.val.photon)
            
        if self.err is None:
            err='No'
        else:
            err=str("%e"%self.err.photon)

        range1=str("%.3f"%self.idx_range[0])
        
        range2=str("%.3f"%self.idx_range[1])
        
        print ("---> name = %-16s range=[%-6s,%-6s] log(Hz)  photon.val=%-13s, err=%-13s"%(self.name,range1,range2,val,err))

class index_array(object):
    """
    Class to handle an array of :class:`index` objects
    """
    def __init__(self):
        self.idx_array=[]
    
    def add_index(self,name=None,data_type=None,val=None,err=None,idx_range=[]):
        self.idx_array.append(index(name=name,data_type=data_type,val=val,err=err,idx_range=idx_range))
    
    def get_by_name(self,name):
        for pi in range(len(self.idx_array)):
            if self.idx_array[pi].name==name:
                return self.idx_array[pi]
        else:
            print ("no index with name %s found"%name)
            return None
            
    def show_pars(self):
        for par in self.idx_array:
            par.show_val()
#----------------------------------------------------

#----------------------------------------------------
def spectral_index_range(name):
        spectral_range_dic={}
        spectral_range_dic['radio']=[6,10]
        spectral_range_dic['radio_mm']=[10,11]
        spectral_range_dic['mm_IR']=[11,13]
        spectral_range_dic['IR_Opt']=[13,14]
        spectral_range_dic['Opt_UV']=[14,16]
        spectral_range_dic['UV_X']=[15,17.5]
        spectral_range_dic['BBB']=[15,16]
        spectral_range_dic['X']=[16,19]
        spectral_range_dic['Fermi']=[22.38,25.38]
        spectral_range_dic['TeV'] = [25.00, 28.38]
        return spectral_range_dic[name]
        
        
def sync_fit_range(name,indices):
        spectral_range_dic={}
        spectral_range_dic['blind']=[9,19]        
        spectral_range_dic['ISP']=[10,17]
        spectral_range_dic['LSP']=[10,17]
        spectral_range_dic['HSP']=[11,20]

        if name=='ISP' and indices.get_by_name('X').val is not None  :
            if  indices.get_by_name('X').val.sed<0:
                return [10,18]
        else:
            return spectral_range_dic[name]
        
      

        return spectral_range_dic[name]
        
def IC_fit_range(name):
    """

    :param name:
    :return:
    """
    spectral_range_dic={}
    spectral_range_dic['blind']=[16,28]
    spectral_range_dic['ISP']=[17,28]
    spectral_range_dic['LSP']=[17,28]
    spectral_range_dic['HSP']=[22,28]

    return spectral_range_dic[name]





class peak_values(object):
        """
        This Class is designed to store the  SED peak values
    
        
        """
        
        def __init__(self,name,nu_p_val=None,nu_p_err=None,nuFnu_p_val=None,nuFnu_p_err=None,curv_val=None,curv_err=None):
            """
            Constructor

            members:

                name->the name member is used to refer to Synchrortorn (S) or Inverse Compoton (IC)
                      components

                the other members are slef-esplicative
            """
            
            self.name=name 
            self.nu_p_err=nu_p_err
            self.nu_p_val=nu_p_val
            self.nuFnu_p_err=nuFnu_p_err
            self.nuFnu_p_val=nuFnu_p_val
            self.curvature=curv_val
            self.curvature_err=curv_err

            
        def update(self,fit_model,fit_law_name):
            #!!CHANGE TO GET PAR BY TYPE
            """ the show() method returns the class members values """
            #print('-->',fit_model,fit_law_name)
            #print(fit_model.show_model())
            Ep=fit_model.parameters.get_par_by_name(fit_law_name,'Ep')
            Sp=fit_model.parameters.get_par_by_name(fit_law_name,'Sp')
            curv=fit_model.parameters.get_par_by_name(fit_law_name,'b')
            #print('-->',Ep,type(Ep),dir(Ep))
            self.nu_p_val=Ep.best_fit_val
            self.nu_p_err=Ep.best_fit_err
            self.nuFnu_p_val=Sp.best_fit_val
            self.nuFnu_p_err=Sp.best_fit_err
            self.curvature=curv.best_fit_val
            self.curvature_err=curv.best_fit_err
            
            
        
        def show(self):
            if self.nu_p_val is None:
                nu_p_val='No'
            else:
                nu_p_val=str("%+e"%self.nu_p_val)
            
            
            if self.nu_p_err is None:
                nu_p_err='No'
            else:
                nu_p_err=str("%+e"%self.nu_p_err)
            
            if self.nuFnu_p_val is None:
                nuFnu_p_val='No'
            else:
                nuFnu_p_val=str("%+e"%self.nuFnu_p_val)
            
            
            if self.nuFnu_p_err is None:
                nuFnu_p_err='No'
            else:
                nuFnu_p_err=str("%+e"%self.nuFnu_p_err)
            
            
            if self.curvature is None:
                curvature='No'
            else:
                curvature=str("%+e"%self.curvature)
            
            
            
            if self.curvature_err is None:
                curvature_err='No'
            else:
                curvature_err=str("%+e"%self.curvature_err)
            
            """ the update() method returns update the class members"""
            
            print ("---> %-10s nu_p=%-13s (err=%-13s)  nuFnu_p=%-13s (err=%-13s) curv.=%-13s (err=%-13s)"%(self.name,nu_p_val,nu_p_err,nuFnu_p_val,nuFnu_p_err,curvature,curvature_err))
        #def update(nu_p_val=None,nu_p_err=None,nuFnu_p_val=None,nuFnu_p_err=None,curv=None,curv_err=None):
            
    
#----------------------------------------------------


#
class SEDShape(object):
    """
    
    This handle the SED shaping process

    """    
    
    def __init__(self,sed_data):
        
        self.indices=index_array()
        self.indices.add_index(name='radio',data_type='sed')    
        self.indices.add_index(name='radio_mm',data_type='sed')   
        self.indices.add_index(name='mm_IR',data_type='sed')
        self.indices.add_index(name='IR_Opt',data_type='sed')
        self.indices.add_index(name='Opt_UV',data_type='sed')
        self.indices.add_index(name='BBB',data_type='sed')
        self.indices.add_index(name='UV_X',data_type='sed')
        self.indices.add_index(name='X',data_type='sed')
        self.indices.add_index(name='Fermi',data_type='sed')
        self.indices.add_index(name='TeV', data_type='sed')
        
        self.sed_data=sed_data
        
        self.cosmo=sed_data.cosmo

       
        self.BBB=None
        self.disk=None
        self.L_host=None
        self.L_host_err=None
        
        self.L_Disk=None
        self.L_Disk_err=None
        
        
        self.T_Disk=None
        
        self.nu_p_Disk=None
        self.nu_p_Disk_err=None
       
        self.S_nu_max=None
        self.IC_nu_max=None
        
        self.obj_class=None
        
        
        

        self.S_peak=peak_values(name='sync')
        self.S_LE_slope=None
        
        self.IC_peak=peak_values(name='IC')
        self.IC_LE_slope=None

        self.host_gal=None
        self.sync_fit_model=None
        self.IC_fit_model=None

    def find_class(self,E_S):
        """

        method to evaluate obj class 'L/I/HPS' according
        to Ep
        
        args: E_S
        
        the obj_class member is updated
        
        obj_class==None means undefined class

        """

        if E_S>6 and E_S<14:
     
            self.obj_class='LSP'
            print ("--> class: ", self.obj_class)
        elif E_S>=14 and  E_S<=15:
             
            self.obj_class='ISP'
            print ("---> class: ", self.obj_class)
        elif E_S>15:
             
            self.obj_class='HSP'
            print ("---> class: ", self.obj_class)
        else:
            self.obj_class=None
            print ("---> undefined class for src",self.obj_class)
    
    
    def show_values(self):
        
        print (section_separator)
        
        print ('*** SEDShape values ***')

        print ('---> spectral inidces values')
        
        for index in self.indices.idx_array:
        
            index.show_val()
        
        print()
        print()
        print ("---> S/IC peak values")
        
        self.S_peak.show()
        
        print ()
        
        self.IC_peak.show()

        print ()
        print ()

        if self.L_Disk is not None  :
            
            print ("---> Disk peak values")
            print ("---> %-10s nu_p_Diks=%+13e (err=%+13e)  L_Disk=%+13e (err=%+13e) T_Disk=%+13e "%('Disk',self.nu_p_Disk,
                                                                                                           self.nu_p_Disk_err,
                                                                                                           self.L_Disk,
                                                                                                           self.L_Disk_err,
                                                                                                           self.T_Disk))
        
        print (section_separator)
    
    
    
    def save_values(self,name):

        _v=[]

        _dt=[('src_name','S32')]
        _v = [name]
        for index in self.indices.idx_array:
            _dt.append((index.name, 'f8'))
            _dt.append((index.name + '_err', 'f8'))
            if index.val is not None:

                _v.extend([index.val.photon, index.err.photon])
            else:
                _v.extend([None,None])



        _dt.append(('nu_p_S', 'f8'))
        _dt.append(('nu_p_S_err', 'f8'))
        _v.extend([self.S_peak.nu_p_val,self.S_peak.nu_p_err])

        _dt.append(('nuFnu_p_S', 'f8'))
        _dt.append(('nuFnu_p_S_err', 'f8'))
        _v.extend([self.S_peak.nuFnu_p_val, self.S_peak.nuFnu_p_err])

        _dt.append(('nu_p_IC', 'f8'))
        _dt.append(('nu_p_IC_err', 'f8'))
        _v.extend([self.IC_peak.nu_p_val, self.IC_peak.nu_p_err])

        _dt.append(('nuFnu_p_IC', 'f8'))
        _dt.append(('nuFnu_p_IC_err', 'f8'))
        _v.extend([self.IC_peak.nuFnu_p_val, self.IC_peak.nuFnu_p_err])


        _out=np.zeros(1,dtype=_dt)
        _out[0]=tuple(_v)

        t=Table(_out)
        t.write(name,overwrite=True)

    def eval_indices(self,minimizer='lsb',silent=True,show_fit_report=False):
        """
        
        This methods evaluates the indices for the SED
        indices are istances of the index_array () class
        
        """
        
        
        print (section_separator)
        
        print ("*** evaluating spectral indices for data ***")
        
        self.index_models=[]
        
        for index in self.indices.idx_array:
            do_fit=self.check_adapt_range_size(self.sed_data.data['nu_data_log'],index,3,silent=silent)
            if do_fit==True:
                loglog_poly=LogLinear()
                loglog_pl=FitModel(cosmo=self.cosmo, name='%s'%index.name,loglog_poly=loglog_poly)
                #print(10.**index.idx_range[0],10.**index.idx_range[1])
                mm,best_fit=fit_SED(loglog_pl,
                                 self.sed_data,
                                 10.**index.idx_range[0],
                                 10.**index.idx_range[1],
                                 loglog=True,
                                 silent=silent,
                                 fitname='spectral-indices-best-fit',
                                 minimizer=minimizer,
                                 use_UL=True)
                #val,err=do_linear_fit(self.SEDdata.nu_data_log,self.SEDdata.nuFnu_data_log,dy=self.SEDdata.dnuFnu_data_log,x_range=index.idx_range)




                par=loglog_pl.parameters.get_par_by_name(loglog_poly.name,'alpha')

                index.assign_val(val = par.best_fit_val, err =par.best_fit_err)

                if silent is False:
                    index.show_val()

                self.index_models.append(loglog_pl)

                if show_fit_report==True and silent is False:
                    best_fit.show_report()
                    print()

            if silent is False:
                print()

        print (section_separator)
            

    def plot_indices(self,plot_obj=None):
        if plot_obj is None:
            plot_obj=PlotSED(sed_data=self.sed_data)

        for model in self.index_models:
            plot_obj.add_model_plot(model,label=model.name,line_style='--')
        
        return plot_obj

    def plot_shape_fit(self, plot_obj=None):
        if plot_obj is None:
            plot_obj = PlotSED(sed_data=self.sed_data)

            if self.sync_fit_model is not None:
                plot_obj.add_model_plot(self.sync_fit_model, label='sync, poly-fit')

            if self.host_gal is not None:
                plot_obj.add_model_plot(self.host_gal, label='host-gal')

            if self.BBB is not None:
                plot_obj.add_model_plot(self.BBB, label='BBB')

            if self.disk is not None:
                plot_obj.add_model_plot(self.disk, label='disk')

            #plot_obj.add_model_plot(self.sync_fit_model, label='sync+host, poly-fit')
            if self.IC_fit_model is not None:
                plot_obj.add_model_plot(self.IC_fit_model, label='IC, poly-fit')

        return plot_obj

    def sync_fit(self,check_host_gal_template=False,
                 check_BBB_template=False,
                 check_disk=False,
                 fit_range=None,
                 nu_min=None,
                 nu_max=None,
                 Ep_start=None,
                 use_log_par=False,
                 minimizer='lsb',
                 silent=True,
                 show_fit_report=False):
        
        """
        This method analyses the synchrotron shape by means
        of log-log polynomial fits

        The following paremeters are estimated:

         1) the SED peak frequency Ep
         2) the SED peak flux Sp
         3) the curvature at the peak
         4) checks for the host galaxy
         
         -) first a log-log cubic fit is performed, with the 'blind' interval

         -) the SED class 'I/L/HPS' is set according to Ep, by find_class

         -) the fit range is changed from 'blind', according to the
            value returned by find_class, using the function sync_fit_range

         -) the SED nu, and nuFnu are generated for the 'blind' fit

         -) if the option is selected in the call of the method
            the estimate of the host galaxy is performed            

         -) a second run improve the obj class is performed

        
        values are stored in the class :class: 'peak_values'

        """
        
        print  (section_separator)
                  
        print ("*** Log-Polynomial fitting of the synchrotron component ***")
        
    
       
        if use_log_par==True:
            fit_law = LogParabolaEp()
        else:
            fit_law=LogCubic()

        self.sync_fit_model=FitModel(cosmo=self.cosmo, loglog_poly=fit_law,name='sync_model')
        self.sync_fit_model.set(fit_law.name,'b',fit_range_max=0)
        self.sync_fit_model.set(fit_law.name,'b',val_max=0)
        self.sync_fit_model.set(fit_law.name,'b',fit_range_min=-10)
        self.sync_fit_model.set(fit_law.name,'b',val_min=-10)

        if Ep_start is not None  :
            self.sync_fit_model.set(fit_law.name,'Ep',val=Ep_start)

        mm, sync_best_fit = self.do_sync_fit(self.sync_fit_model,
                                             fit_law.name,
                                             fit_range=fit_range,
                                             check_disk=check_disk,
                                             check_BBB=check_BBB_template,
                                             check_host=check_host_gal_template,
                                             # use_log_par=False,
                                             Ep_start=Ep_start,
                                             minimizer=minimizer,
                                             silent=silent,
                                             show_fit_report=show_fit_report)
        #print "bbb", self.sync_fit_model.SED.nu

        #Ep=self.S_peak.nu_p_val

        #logparpl=LogParabolaPL()
        #fit_model=FitModel(loglog_poly=logparpl,name='sync_logpar_pl_model')
        #self.do_sync_fit(fit_model,fit_range=fit_range,
        #                 check_disk=check_disk,
        #                 check_BBB=check_BBB_template,
        #                 check_host=check_host_gal_template,
        #                 use_log_par=True,
        #                 Ep=Ep)



        #print "bbb", self.sync_fit_model.SED.nu

        self.sync_best_fit=sync_best_fit

        return mm,sync_best_fit
       
    
    def save_sync_fit_report(self,name=None):
        self.sync_best_fit.save_report(name=name)

    
    def do_sync_fit(self,
                    fit_model,
                    fit_law_name,
                    fit_range=None,
                    check_disk=False,
                    check_BBB=False,
                    check_host=False,
                    #use_log_par=False,
                    Ep_start=None,
                    no_check=False,
                    minimizer='lsb',
                    silent=True,
                    show_fit_report=True):
        
       
        
      
        #fit_model1=copy.deepcopy(fit_model)

        if fit_range is None:
            s_fit_range=sync_fit_range('blind',self.indices)
        else:
            s_fit_range=fit_range
            
        print ("---> first blind fit run,  fit range:",s_fit_range)

        if silent == False:
            fit_model.show_pars()

        if Ep_start is None:
            Ep=find_max_cubic(self.sed_data.data['nu_data_log'],self.sed_data.data['nuFnu_data_log'] ,x_range=s_fit_range)

        else:
            Ep=Ep_start


        #if use_log_par == False:
        fit_model.set(fit_law_name,'Ep', val=Ep)

        #else:
        #    fit_model1.set('E0', val=Ep - 1)

        mm,best_fit=fit_SED(fit_model,self.sed_data,10.**s_fit_range[0],10.**s_fit_range[1],loglog=True,silent=silent,fitname='sync-shape-fit',minimizer=minimizer,use_UL=True)

        self.S_peak.update(fit_model,fit_law_name)
        self.find_class(self.S_peak.nu_p_val)

        refit=False
        if check_disk==True : 
            self.add_disk(fit_model)
            refit=True
        if check_BBB==True:
            self.add_BBB_template(fit_model)
            refit=True
        if check_host==True:
            self.add_host_template(fit_model)
            refit=True

        if silent == False:
            if show_fit_report == True:
                best_fit.show_report()
        print()
        #Ep=None



        if refit==True:
            if fit_range is None:
                s_fit_range = sync_fit_range(self.obj_class, self.indices)

            mm,best_fit=fit_SED(fit_model,self.sed_data,10.**s_fit_range[0],10.**s_fit_range[1],loglog=True,silent=True,fitname='sync-shape-fit',minimizer=minimizer,use_UL=True)

            if silent == False:
                best_fit.show_report()

            #if use_log_par==False:
            self.S_peak.update(fit_model,fit_law_name)

            self.find_class(self.S_peak.nu_p_val)



        print()
        print()
        fit_model.show_best_fit_pars()
        self.S_peak.show()

        self.S_nu_max= self.get_nu_max(self.sed_data.data['nu_data_log'] ,s_fit_range)
        
        
        self.set_S_LE_slope(fit_model,fit_law_name,use_log_par=False)


        if check_disk==True :
            self.disk.set_disk_pars(fit_model,'Disk')
            self.L_Disk=self.disk.L_Disk
            self.L_Disk_err=self.disk.L_Disk_err
            self.T_Disk=self.disk.T_Disk
            self.T_Disk_err=self.disk.T_Disk_err
            self.nu_p_Disk=self.disk.nu_p_Disk
            self.nu_p_Disk_err=self.disk.nu_p_Disk_err


        if check_host==True :
            self.host_gal.set_host_pars(fit_model,'host_galaxy')

        if check_BBB==True :
            self.BBB.set_BBB_pars(fit_model,'BBB')
            self.L_Disk=self.BBB.nuLnu_p_BBB
            self.L_Disk_err=self.BBB.nuLnu_p_BBB_err
            self.T_Disk=self.BBB.T_Disk
            self.T_Disk_err=self.BBB.T_Disk_err
            self.nu_p_Disk=self.BBB.nu_p
            self.nu_p_Disk_err=self.BBB.nu_p_err

        print(section_separator)

        return mm,best_fit
        
    
    def add_disk(self,fit_model):
        
        disk=Disk(cosmo=fit_model.cosmo, z=self.sed_data.z,name='Disk')
        
        fit_model.add_component(disk)

        fit_model.set('Disk','nuFnu_p',val=(10**self.S_peak.nuFnu_p_val),fit_range=[10**(self.S_peak.nuFnu_p_val-2),10**(self.S_peak.nuFnu_p_val+2 )])
        
        fit_model.set('Disk','T_Disk',val=1E4,fit_range=[5000,1E5])
 
        self.disk=disk
 
 
 
 
    
    def   add_BBB_template(self,fit_model):
        
        BBB_template=SpectralTemplateLogLog.template_factory('BBB', cosmo=fit_model.cosmo, z=self.sed_data.z, name='BBB')
        
        fit_model.add_component(BBB_template)
   
        fit_model.set('BBB','nuFnu_p_BBB',val=(self.S_peak.nuFnu_p_val),fit_range=[self.S_peak.nuFnu_p_val-2,self.S_peak.nuFnu_p_val+2 ])
        
        fit_model.set('BBB','nu_scale',val=0,fit_range=[-0.5,0.5])
   
        self.BBB=BBB_template
    
    
    
        
    
    
        

    def add_host_template(self,fit_model):
        
        host_gal=SpectralTemplateLogLog.template_factory('host_galaxy', cosmo=fit_model.cosmo, z=self.sed_data.z, name='host_galaxy')
                 
        fit_model.add_component(host_gal)

        #print('nuFnu_p_host',self.S_peak.nuFnu_p_val)
        fit_model.set('host_galaxy','nuFnu_p_host',val=(self.S_peak.nuFnu_p_val),fit_range=[self.S_peak.nuFnu_p_val-2,self.S_peak.nuFnu_p_val+2 ])
        
        fit_model.set('host_galaxy','nu_scale',val=0,fit_range=[-0.5,0.5])
   
        self.host_gal=host_gal
    
    
    
    
    def set_S_LE_slope(self,fit_func,fit_law_name,use_log_par):
        if use_log_par==False:
            self.S_LE_slope=2.0
        else:
            self.S_LE_slope=fit_func.parameters.get_par_by_name(fit_law_name,'alpha').best_fit_val

    
    
    
    
      
   
    
    def IC_fit(self,fit_range=None,use_log_par=False,Ep_start=None,minimizer='minuit',silent=False):
    
        print  (section_separator)
              
        print ("*** Log-Polynomial fitting of the IC component ***")
        
        if fit_range is None:
            fit_range=IC_fit_range(self.obj_class)
        
        print ("---> fit range:",fit_range)
        

        


        if Ep_start is None:
            Ep=find_max_cubic(self.sed_data.data['nu_data_log'],self.sed_data.data['nuFnu_data_log'] ,x_range=fit_range)
        else:
            Ep=Ep_start


        if use_log_par==True:
            print("---> LogParabola fit")
            logpar_model = LogParabolaEp()
            fit_law_name = logpar_model.name
            self.IC_fit_model = FitModel(cosmo=self.cosmo, loglog_poly=logpar_model, name='IC_log-par_model')



            if Ep is not None:
                self.IC_fit_model.set(fit_law_name,'Ep', val=Ep)

            if silent == False:
                self.IC_fit_model.show_pars()

            mm,best_fit = fit_SED(self.IC_fit_model, self.sed_data, 10. ** fit_range[0], 10. ** fit_range[1], loglog=True,
                               silent=True, fitname='IC-shape-fit', minimizer=minimizer,use_UL=True)

            if silent == False:
                best_fit.show_report()

            best_fit_model = self.IC_fit_model
        else:
            print("---> LogCubic fit")
            cubic_model = LogCubic()
            fit_law_name=cubic_model.name
            self.IC_fit_model = FitModel(cosmo=self.cosmo, loglog_poly=cubic_model, name='IC_log-cubic_model')

            self.IC_fit_model.set(fit_law_name,'b', fit_range_max=0)
            self.IC_fit_model.set(fit_law_name,'b', val_max=0)
            self.IC_fit_model.set(fit_law_name,'b', fit_range_min=-10)
            self.IC_fit_model.set(fit_law_name,'b', val_min=-10)
            if Ep is not None:
                self.IC_fit_model.set(fit_law_name,'Ep', val=Ep)

            try:
                mm,best_fit=fit_SED(self.IC_fit_model,self.sed_data,10.**fit_range[0],10.**fit_range[1],loglog=True,silent=True,fitname='IC-shape-fit',minimizer=minimizer,use_UL=True)

                if silent == False:
                    best_fit.show_report()

                best_fit_model=self.IC_fit_model


            except Exception as e:

                print ("---> LogCubic fit failed",e)
                print ("---> try LogParabola ")
                logpar_model=LogParabolaEp()
                fit_law_name=logpar_model.name
                self.IC_fit_model=FitModel(cosmo=self.cosmo, loglog_poly=logpar_model,name='IC_log-par_model')
                if Ep is not None:
                    self.IC_fit_model.set(fit_law_name,'Ep', val=Ep)

                if silent == False:
                    self.IC_fit_model.show_pars()

                mm,best_fit=fit_SED(self.IC_fit_model,self.sed_data,10.**fit_range[0],10.**fit_range[1],loglog=True,silent=True,fitname='IC-shape-fit',minimizer=minimizer,use_UL=True)

                if silent == False:
                    best_fit.show_report()

                best_fit_model=self.IC_fit_model
            
            
            
        self.IC=best_fit_model
        
        self.IC_peak.update(best_fit_model,fit_law_name)

        print()
        print()
        self.IC_fit_model.show_best_fit_pars()

        self.IC_peak.show()
        
      
        self.IC_nu_max= self.get_nu_max(self.sed_data.data['nu_data_log'],fit_range)
        
        print  (section_separator)

        return mm,best_fit
    
    
    def get_nu_max(self,nu,fit_range):

        msk = filter_interval(nu,fit_range)
        return nu[msk].max()
        
        
    def check_adapt_range_size(self,x,index,min_size,silent=False):
        do_fit=True
        x_range=[index.idx_range[0],index.idx_range[1]]
        msk = filter_interval(x,x_range)
        x1=x[msk]
        delta=0.1
        delta_tot=0
        if silent is False:
            print("---> initial range for index %s  set to [%f,%f]" % (index.name, index.idx_range[0], index.idx_range[1]))
        while len(x1)<min_size and delta_tot<1.0:
            delta_tot+=delta
            x_range[0]-=delta
            x_range[1]+=delta
            msk = filter_interval(x,x_range)
            x1=x[msk]

        if len(x1)>=min_size:
            do_fit=True
            index.idx_range=x_range
            if silent is False:
                print ("---> range for index %s updated  to [%f,%f]"%(index.name,index.idx_range[0],index.idx_range[1] ))
        else:
            do_fit=False
            if silent is False:
                print("---> not enough data in range for index%s " % (index.name))

        return  do_fit
    
def find_E0(b,a,Ep):
    """returns the value of E0  for
    a log_par+pl distribution
    
    Args:
        b:  curvature
        a:  spectral index in the PL branch
        Ep: peak value 

    Returns:
        
        
    """
    if (b!=0.0):
        c=(3-a)/(-2*b)
        print ("Ep=",Ep)
        return Ep/(pow(10,c))
    else:
        return Ep

#----------------------------------------------------   