__author__ = "Andrea Tramacere"


#Peak values
Sync_nuFnu_p_dic={'obs':'Sync.spec.nuFnu_peak_obs'}
Sync_nuFnu_p_dic['src']='Sync.spec.nuLnu_peak_src'
Sync_nuFnu_p_dic['blob']='Sync.spec.nuLnu_peak_blob'


Sync_nu_p_dic={'obs':'Sync.spec.nu_peak_obs'}
Sync_nu_p_dic['src']='Sync.spec.nu_peak_src'
Sync_nu_p_dic['blob']='Sync.spec.nu_peak_blob'


SSC_nuFnu_p_dic={'obs':'SSC.spec.nuFnu_peak_obs'}
SSC_nuFnu_p_dic['src']='SSC.spec.nuLnu_peak_src'
SSC_nuFnu_p_dic['blob']='SSC.spec.nuLnu_peak_blob'


SSC_nu_p_dic={'obs':'SSC.spec.nu_peak_obs'}
SSC_nu_p_dic['src']='SSC.spec.nu_peak_src'
SSC_nu_p_dic['blob']='SSC.spec.nu_peak_blob'

# nu_src_start_stop_dict={'Sync':['nu_start_Sync', 'nu_start_Sync']}
# nu_src_start_stop_dict['SSC']=['nu_start_SSC', 'nu_stop_SSC']
# nu_src_start_stop_dict['EC_BLR']=['', '']
# nu_src_start_stop_dict['EC_DT']=['', '']
# nu_src_start_stop_dict['EC_Disk']=['', '']
# nu_src_start_stop_dict['EC_CMB']=['', '']
# nu_src_start_stop_dict['EC_CMB_stat']=['', '']
# nu_src_start_stop_dict['Bremss_ep']=['', '']
# nu_src_start_stop_dict['PP_gamma']=['', '']
# nu_src_start_stop_dict['PP_neutrino_tot']=['', '']
# nu_src_start_stop_dict['PP_neutrino_e']=['', '']
# nu_src_start_stop_dict['PP_neutrino_mu']=['', '']






#Spectral components
nuFnu_obs_dict={'Sum':['core.nuFnu_sum_grid', 'core.nu_grid']}
nuFnu_obs_dict['Sync']=['Sync.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['SSC']=['SSC.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_BLR']=['BLR.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['BLR']=['BLR.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_DT']=['DT.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['DT']=['DT.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_Corona']=['Corona.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['Corona']=['Corona.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['Star']=['Star.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_Star']=['Star.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_Disk']=['Disk.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['Disk']=['Disk.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_CMB']=['CMB.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['EC_CMB_stat']=['CMB.ec.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['Bremss_ep']=['Bremss_ep.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['PP_gamma']=['PP_gamma.spec.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['PP_neutrino_tot']=['PP_neutrino.spec_tot.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['PP_neutrino_e']=['PP_neutrino.spec_e.nuFnu_grid', 'core.nu_grid']
nuFnu_obs_dict['PP_neutrino_mu']=['PP_neutrino.spec_mu.nuFnu_grid', 'core.nu_grid']






#seed-photon fields components
n_seed_dic={'DT':['DT.spec.n_nu','DT.spec.nu']}
n_seed_dic['EC_DT']=['DT.spec.n_nu','DT.spec.nu']
n_seed_dic['Corona']=['Corona.spec.n_nu','Corona.spec.nu']
n_seed_dic['EC_Corona']=['Corona.spec.n_nu','Corona.spec.nu']

n_seed_dic['EC_Disk']=['Disk.spec.n_nu','Disk.spec.nu']
n_seed_dic['Disk']=['Disk.spec.n_nu','Disk.spec.nu']
n_seed_dic['EC_BLR']=['BLR.spec.n_nu_DRF','BLR.spec.nu_DRF']
n_seed_dic['EC_CMB']=['CMB.spec.n_nu','CMB.spec.nu']
n_seed_dic['EC_S']=['CMB.spec.n_nu','CMB.spec.nu']
n_seed_dic['SSC']=['Sync.spec.n_nu','Sync.spec.nu']
n_seed_dic['Star']=['Star.spec.n_nu','Sync.spec.nu']


#nuLnu_dic={'SUM':['','nu_']}
#nuLnu_dic['Sync']=['nuFnu_Sync','nu_']
#nuLnu_dic['SSC']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_BLR']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_DT']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_Disk']=['nuFnu_Sync','nu_']





#Electron distributions

gamma_dic_e={'electron_distr':['emitters.Ne_jetset','emitters.griglia_gamma_Ne_log']}
gamma_dic_e_equilibrium={'e_inj':['emitters.Q_inj_e_primaries','emitters.griglia_gamma_Ne_log']}

gamma_dic_p={'proton_distr':['emitters.Np_jetset','emitters.griglia_gamma_Np_log']}
gamma_dic_pp_e_second={'e_second_inj':['emitters.Q_inj_e_second','emitters.griglia_gamma_Ne_log']}

s_dic={'pl':'p'}
s_dic['lppl']='s'
s_dic['lp']='s'
s_dic['bkn']='p'
s_dic['plc']='p'

s1_dic={'bkn':'p_1'}

gamma_cut_dic={'pl':'gmax'}
gamma_cut_dic['lppl']='gamma0_log_parab'
gamma_cut_dic['lp']='gmax'
gamma_cut_dic['plc']='gamma_cut'

r_dic={'pl':'r'}
r_dic['lppl']='r'
r_dic['lp']='r'
r_dic['lpep']='r'

gamma_3p_dic={'lpep':'gammap_log_parab'}
gamma_3p_dic['bkn']='gamma_break'


available_N_distr=['lp', 'pl', 'lppl', 'lpep', 'plc', 'bkn', 'spitkov', 'lppl_pile_up', 'bkn_pile_up']
available_N_distr_descr=['log-parabola',
                         'powerlaw',
                         'log-parabola with low-energy powerlaw branch',
                         'log-parabola defined by peak energy',
                         'powerlaw with cut-off',
                         'broken powerlaw',
                         'spitkov',
                         'log-parabola with low-energy powerlaw branch and pile-up',
                         'broken powerlaw and pileup']


available_emitters_type=['electrons','protons']#,'electrons-equilibrium']

N_distr_descr={}
for m,d in zip(available_N_distr,available_N_distr_descr):
    N_distr_descr[m]=d


allowed_disk_type=['BB','MultiBB','Mono']
