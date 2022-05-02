import  pytest
from .base_class import TestBase
import numpy as np

class TestTempEv(TestBase):

    def integration_suite(self,plot=False):
        self._all(plot=plot)

    def test_emp_ev_two_zone_cooling_acc(self, plot=False):
        print('--------> test_emp_ev_two_zone_cooling_acc',plot)
        from jetset.jet_emitters_factory import InjEmittersFactory
        from jetset.jet_model import Jet
        jet_model=Jet()
        q_inj=InjEmittersFactory().create_inj_emitters('pl',emitters_type='electrons',normalize=True)
        q_inj.parameters.gmin.val=9
        q_inj.parameters.gmax.val=10
        q_inj.parameters.p.val=0.5

        jet_model.parameters.beam_obj.val=30
        jet_model.parameters.B.val=0.2
        jet_model.parameters.z_cosm.val=0.03
        jet_model.parameters.R.val=5E15

        flare_duration=1.0E5
        duration=flare_duration*10
        t_D0=1.5E5
        t_A0=2.5E4
        T_esc_rad=1E60
        L_inj=5.0E39
        E_acc_max=4E60
        Delta_R_acc_ratio=0.1
        B_ratio=1.0
        T_SIZE=2E4
        NUM_SET=500
        Diff_Index=2.0
        Acc_Index=1.0

        from jetset.jet_timedep import JetTimeEvol
        temp_ev_acc=JetTimeEvol(jet_rad=jet_model,Q_inj=q_inj,inplace=True)

        temp_ev_acc.rad_region.jet.nu_min=1E8
        temp_ev_acc.acc_region.jet.nu_min=1E8
        T_SIZE=np.int(T_SIZE)

        if Delta_R_acc_ratio is not None:
            temp_ev_acc.parameters.Delta_R_acc.val=temp_ev_acc.parameters.R_rad_start.val*Delta_R_acc_ratio

        T_esc_acc=t_A0/(temp_ev_acc.parameters.Delta_R_acc.val/3E10)*2



        temp_ev_acc.parameters.duration.val=duration
        temp_ev_acc.parameters.TStart_Acc.val=0
        temp_ev_acc.parameters.TStop_Acc.val=flare_duration
        temp_ev_acc.parameters.TStart_Inj.val=0
        temp_ev_acc.parameters.TStop_Inj.val=flare_duration
        temp_ev_acc.parameters.T_esc_acc.val=T_esc_acc
        temp_ev_acc.parameters.T_esc_rad.val=T_esc_rad
        temp_ev_acc.parameters.t_D0.val=t_D0
        temp_ev_acc.parameters.t_A0.val=t_A0
        temp_ev_acc.parameters.Esc_Index_acc.val=Diff_Index-2
        temp_ev_acc.parameters.Esc_Index_rad.val=0
        temp_ev_acc.parameters.Acc_Index.val=Acc_Index
        temp_ev_acc.parameters.Diff_Index.val=Diff_Index
        temp_ev_acc.parameters.t_size.val=T_SIZE
        temp_ev_acc.parameters.num_samples.val=NUM_SET
        temp_ev_acc.parameters.E_acc_max.val=E_acc_max
        temp_ev_acc.parameters.L_inj.val=L_inj


        temp_ev_acc.parameters.gmin_grid.val=1.0
        temp_ev_acc.parameters.gmax_grid.val=1E8
        temp_ev_acc.parameters.gamma_grid_size.val=1500

        temp_ev_acc.parameters.B_acc.val=temp_ev_acc.rad_region.jet.parameters.B.val*B_ratio
        temp_ev_acc.init_TempEv()


        only_injection=True
        do_injection=True
        plot_fit_model=True
        plot_fit_distr=True
        plot_emitters=True
        plot_lcs=True
        delta_t_out=1000
        eval_cross_time=False
        rest_frame='obs'
        temp_ev_acc.run(only_injection=only_injection,
                        do_injection=do_injection,
                        cache_SEDs_acc=True, 
                        cache_SEDs_rad=True)

        if plot is True:
            p=temp_ev_acc.plot_tempev_emitters(region='rad',loglog=False,energy_unit='gamma',pow=0)
            p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_A, ls='--')
            p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_DA, ls='--')
            p.setlim(x_max=1E7,x_min=1,y_min=1E-18,y_max=100)

            p=temp_ev_acc.plot_tempev_emitters(region='acc',loglog=False,energy_unit='gamma',pow=0)
            p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_A, ls='--')
            p.ax.axvline(temp_ev_acc.temp_ev.gamma_eq_t_DA, ls='--')
            p.setlim(x_max=1E7,x_min=1,y_min=1E-30,y_max=100)

            p=temp_ev_acc.plot_tempev_model(region='rad',sed_data=None, use_cached = True)
            p.setlim(y_min=1E-18,x_min=1E7)


        lg=temp_ev_acc.rad_region.make_lc(nu1=2.4E22,nu2=7.2E25,name='gamma',eval_cross_time=False,delta_t_out=100,use_cached=True,frame='obs')

    
        if plot is True:
            import matplotlib.pyplot as plt
            plt.plot(lg['time'],lg['flux'])
            plt.xlabel('time (%s)'%lg['time'].unit)
            plt.ylabel('flux (%s)'%lg['flux'].unit)
        
        temp_ev_acc.save_model('two_zone_rad_acc.pkl')
        temp_ev_acc_1=JetTimeEvol.load_model('two_zone_rad_acc.pkl')

        if plot is True:
            p=temp_ev_acc_1.plot_tempev_model(region='rad',sed_data=None, use_cached = True)
        
        lx=temp_ev_acc_1.rad_region.make_lc(nu1=1E17,nu2=1E18,name='X',eval_cross_time=False,delta_t_out=100,use_cached=True,frame='obs')

        if plot is True:
            plt.plot(lx['time'],lx['flux'])
            plt.xlabel('time (%s)'%lg['time'].unit)
            plt.ylabel('flux (%s)'%lg['flux'].unit)

        print('-------->',plot)