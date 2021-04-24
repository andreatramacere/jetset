__author__ = "Andrea Tramacere"


#Peak values
Sync_nuFnu_p_dic={'obs':'nuFnu_peak_Sync_obs'}
Sync_nuFnu_p_dic['src']='nuLnu_peak_Sync_src'
Sync_nuFnu_p_dic['blob']='nuLnu_peak_Sync_blob'


Sync_nu_p_dic={'obs':'nu_peak_Sync_obs'}
Sync_nu_p_dic['src']='nu_peak_Sync_src'
Sync_nu_p_dic['blob']='nu_peak_Sync_blob'


SSC_nuFnu_p_dic={'obs':'nuFnu_peak_SSC_obs'}
SSC_nuFnu_p_dic['src']='nuLnu_peak_SSC_src'
SSC_nuFnu_p_dic['blob']='nuLnu_peak_SSC_blob'


SSC_nu_p_dic={'obs':'nu_peak_SSC_obs'}
SSC_nu_p_dic['src']='nu_peak_SSC_src'
SSC_nu_p_dic['blob']='nu_peak_SSC_blob'

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
nuFnu_obs_dict={'Sum':['nuFnu_sum_grid', 'nu_grid']}
nuFnu_obs_dict['Sync']=['nuFnu_Sync_grid', 'nu_grid']
nuFnu_obs_dict['SSC']=['nuFnu_SSC_grid', 'nu_grid']
nuFnu_obs_dict['EC_BLR']=['nuFnu_EC_BLR_grid', 'nu_grid']
nuFnu_obs_dict['EC_DT']=['nuFnu_EC_DT_grid', 'nu_grid']
nuFnu_obs_dict['DT']=['nuFnu_DT_grid', 'nu_grid']
nuFnu_obs_dict['Star']=['nuFnu_Star_grid', 'nu_grid']
nuFnu_obs_dict['EC_Disk']=['nuFnu_EC_Disk_grid', 'nu_grid']
nuFnu_obs_dict['Disk']=['nuFnu_Disk_grid', 'nu_grid']
nuFnu_obs_dict['EC_CMB']=['nuFnu_EC_CMB_grid', 'nu_grid']
nuFnu_obs_dict['EC_CMB_stat']=['nuFnu_EC_CMB_stat_grid', 'nu_grid']
nuFnu_obs_dict['Bremss_ep']=['nuFnu_bremss_ep_grid', 'nu_grid']
nuFnu_obs_dict['PP_gamma']=['nuFnu_pp_gamma_grid', 'nu_grid']
nuFnu_obs_dict['PP_neutrino_tot']=['nuFnu_pp_neutrino_tot_grid', 'nu_grid']
nuFnu_obs_dict['PP_neutrino_e']=['nuFnu_pp_neutrino_e_grid', 'nu_grid']
nuFnu_obs_dict['PP_neutrino_mu']=['nuFnu_pp_neutrino_mu_grid', 'nu_grid']






#seed-photon fields components
n_seed_dic={'DT':['n_DT','nu_DT']}
n_seed_dic['EC_DT']=['n_DT','nu_DT']

n_seed_dic['EC_Disk']=['n_Disk','nu_Disk']
n_seed_dic['Disk']=['n_Disk','nu_Disk']
n_seed_dic['EC_BLR']=['n_BLR','nu_BLR']
n_seed_dic['EC_CMB']=['n_CMB','nu_CMB']
n_seed_dic['EC_S']=['n_CMB','nu_CMB']
n_seed_dic['SSC']=['n_Sync','nu_Sync']
n_seed_dic['Star']=['n_Star','nu_Sync']


#nuLnu_dic={'SUM':['','nu_']}
#nuLnu_dic['Sync']=['nuFnu_Sync','nu_']
#nuLnu_dic['SSC']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_BLR']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_DT']=['nuFnu_Sync','nu_']
#nuLnu_dic['EC_Disk']=['nuFnu_Sync','nu_']





#Electron distributions

gamma_dic_e={'electron_distr':['Ne','griglia_gamma_Ne_log']}

gamma_dic_p={'proton_distr':['Np','griglia_gamma_Np_log']}
gamma_dic_pp_e_second={'e_second_inj':['Q_inj_e_second','griglia_gamma_Ne_log']}

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


available_emitters_type=['electrons','protons']

N_distr_descr={}
for m,d in zip(available_N_distr,available_N_distr_descr):
    N_distr_descr[m]=d


allowed_disk_type=['BB','MultiBB','Mono']