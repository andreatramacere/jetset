//=========================================================================================
// FILE TEMPLATES
//
//=========================================================================================

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
/**
 * \file funzioni_compton.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief funzioni per il Compton
 *
 */

void flux_header(FILE *fp){
    fprintf(fp, "#######################################################################################################\n");
    fprintf(fp, "#_obs is as observed from the earth: beaming+cosmo   \n");
    fprintf(fp, "#_src is in the AGN rest frame     :cosmo         \n");
    fprintf(fp, "#_blob is in the blob rest frame                   \n");
    fprintf(fp, "# log(nu_obs)|  log(nu_obs*F_nu_obs)|  nu_obs| nu_obs*F_nu_obs|   nu_src| nu_src*L_nu_src|  nu_blob| nu_blob*L_nu_blob|     \n");
    fprintf(fp, "#######################################################################################################\n");
    
}

void flux_DISK_header(FILE *fp){
    fprintf(fp, "#####################################################################################\n");
    fprintf(fp, "#_obs is as observed from the earth: beaming+cosmo               \n");
    fprintf(fp, "#_src is in the AGN rest frame:cosmo                 \n");
    fprintf(fp, "# log(nu_obs)|  log(nu_obs*F_nu_obs)|  nu_obs| nu_obs*F_nu_obs|   nu_src| nu_src*L_nu_src|      \n");
    fprintf(fp, "####################################################################################\n");
    
}

void distr_e_header(FILE *fp){
    fprintf(fp, "#################################################################\n");
    fprintf(fp, "# log(gamma)     log(N(gamma))   gamma     N(gamma)  TeV N(TeV) # \n");
    fprintf(fp, "#################################################################\n");
}

void somma_header(FILE *fp){
    fprintf(fp, "########################################################################\n");
    fprintf(fp, "#_obs is as observed from the earth: beaming+cosmo   \n");
    fprintf(fp, "#nu_obs  nuFnu_somma_obs  nuFnu_Sync_obs nuFnu_comp_obs nu_Fnu_Disk_obs  nu_Fnu_DT_obs nu_Fnu_Star_obs nu_Fnu_EC_Disk_obs nu_Fnu_EC_BLR_obs nu_Fnu_EC_DT_obs  nu_Fnu_EC_Star_obs nu_Fnu_EC_CMB_obs\n");
    fprintf(fp, "########################################################################\n");
}

void somma_log_log_header(FILE *fp){
    fprintf(fp, "########################################################################\n");
    fprintf(fp, "#_obs is as observed from the earth: beaming+cosmo   \n");
    fprintf(fp, "#log(nu_obs) log(nuFnu_somma_obs)  log(nuFnu_Sync_obs) log(nuFnu_SSC_obs)  log(nu_Fnu_Disk_obs) log(nu_Fnu_DT_obs)  log(nu_Fnu_Star_obs) log(nu_Fnu_EC_Disk_obs) log(nu_Fnu_EC_BLR_obs) log(nu_Fnu_EC_DT_obs)   log(nu_Fnu_EC_Star_obs) log(nu_Fnu_EC_CMB_obs)  \n");
    fprintf(fp, "########################################################################\n");
}

void somma_log_log_src_header(FILE *fp){
    fprintf(fp, "########################################################################\n");
    fprintf(fp, "#_src is in the AGN rest frame:cosmo                 \n");
    fprintf(fp, "#nu_src  nuFnu_somma_src  nuFnu_Sync_src nuFnu_SSC_src nu_Fnu_Disk_src  nu_Fnu_DT_src nu_Fnu_Star_src nu_Fnu_EC_Disk_src nu_Fnu_EC_BLR_src nu_Fnu_EC_DT_src nu_Fnu_EC_Star_src nu_Fnu_EC_CMB_src  \n");
    fprintf(fp, "########################################################################\n");
}


