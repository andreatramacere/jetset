#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"


/*
void build_photons(struct spettro *pt_base){

// Grid
    alloc_photons(&(pt_base->core.nu_grid),pt_base->core.nu_grid_size );
    alloc_photons(&(pt_base->core.nuFnu_sum_grid),pt_base->core.nu_grid_size);
    //pt_base->core.nu_grid = malloc(200 * sizeof (double));
    //pt_base->core.nuFnu_sum_grid = malloc(200 * sizeof (double));
    //double nu_grid[static_spec_arr_size];
    //double nuFnu_sum_grid[static_spec_arr_size];





//Synch
//  double nuFnu_Sync_grid[static_spec_arr_size];
//  double j_Sync[static_spec_arr_size];
//  double alfa_Sync[static_spec_arr_size];
//  double I_nu_Sync[static_spec_arr_size];
//  double nu_Sync[static_spec_arr_size];
//  double nu_Sync_obs[static_spec_arr_size];
//  double n_Sync[static_spec_arr_size];
//  double nuF_nu_Sync_obs[static_spec_arr_size];

    alloc_photons(&(pt_base->Sync.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->Sync.spec.j_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.alfa_Sync),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.spec.I_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.spec.nu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.spec.n_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Sync.spec.nuFnu_obs),pt_base->core.nu_seed_size);



//Hadrons
//  double j_pp[static_spec_arr_size];
//  double nu_pp[static_spec_arr_size];
//  double nuF_nu_pp_obs[static_spec_arr_size];
//
    alloc_photons(&(pt_base->j_pp),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->nu_pp),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->nuF_nu_pp_obs),pt_base->core.nu_seed_size);



//SSC
//  double nuFnu_SSC_grid[static_spec_arr_size];
//  double q_comp[static_spec_arr_size];
//  double j_comp[static_spec_arr_size];
//  double j_EC[static_spec_arr_size];
//  double nu_SSC[static_spec_arr_size];
//  double nu_SSC_obs[static_spec_arr_size];
//  double nuF_nu_SSC_obs[static_spec_arr_size];

    alloc_photons(&(pt_base->SSC.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->SSC.q_comp),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->SSC.spec.j_nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->j_EC),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->SSC.spec.nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->SSC.spec.nu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->SSC.spec.nuFnu_obs),pt_base->core.nu_IC_size);


//EC star
//  double nuFnu_EC_Star_grid[static_spec_arr_size];
//	double nu_EC_Star[static_spec_arr_size];
//	double nu_EC_Star_obs[static_spec_arr_size];
//	double nuF_nu_EC_Star_obs[static_spec_arr_size];



    alloc_photons(&(pt_base->Star.ec.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->Star.ec.spec.nuFnu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->Star.ec.spec.nu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->Star.ec.spec.nu),pt_base->core.nu_IC_size);

//  double nuFnu_Star_grid[static_spec_arr_size];
//  double I_nu_Star[static_spec_arr_size];
//	double J_nu_Star_disk_RF[static_spec_arr_size];
//	double I_nu_Star_disk_RF[static_spec_arr_size];
//	double nu_Star[static_spec_arr_size];
//	double nu_Star_obs[static_spec_arr_size];
//	double nu_Star_disk_RF[static_spec_arr_size];
//	double nuF_nu_Star_obs[static_spec_arr_size];
//	double n_Star[static_spec_arr_size];


    alloc_photons(&(pt_base->Star.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->Star.spec.I_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.J_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.I_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.nu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.nuFnu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Star.spec.n_nu),pt_base->core.nu_seed_size);







//EC CMB
//  double nuFnu_EC_CMB_grid[static_spec_arr_size];
//	double nu_EC_CMB[static_spec_arr_size];
//	double nu_EC_CMB_obs[static_spec_arr_size];
//	double nuF_nu_EC_CMB_obs[static_spec_arr_size];

    alloc_photons(&(pt_base->CMB.ec.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->CMB.ec.spec.nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->CMB.ec.spec.nu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->CMB.ec.spec.nuFnu_obs),pt_base->core.nu_IC_size);

//	double I_nu_CMB[static_spec_arr_size];
//	double I_nu_CMB_disk_RF[static_spec_arr_size];
//	double nu_CMB[static_spec_arr_size];
//	double nu_CMB_disk_RF[static_spec_arr_size];
//	double n_CMB[static_spec_arr_size];

    alloc_photons(&(pt_base->CMB.spec.I_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->CMB.spec.I_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->CMB.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->CMB.spec.nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->CMB.spec.n_nu),pt_base->core.nu_seed_size);



//
//
//EC CMB stat
//  double nuFnu_EC_CMB_stat_grid[static_spec_arr_size];
//  double nu_EC_CMB_stat[static_spec_arr_size];
//  double nu_EC_CMB_stat_obs[static_spec_arr_size];
//  double nuF_nu_EC_CMB_stat_obs[static_spec_arr_size];

//    alloc_photons(&(pt_base->nuFnu_EC_CMB_stat_grid),pt_base->core.nu_grid_size);

//    alloc_photons(&(pt_base->nuF_nu_EC_CMB_stat_obs),pt_base->core.nu_IC_size);
//    alloc_photons(&(pt_base->nu_EC_CMB_stat),pt_base->core.nu_IC_size);
//    alloc_photons(&(pt_base->nu_EC_CMB_stat_obs),pt_base->core.nu_IC_size);

//  double I_nu_CMB_stat[static_spec_arr_size];
//  double I_nu_CMB_disk_RF_stat[static_spec_arr_size];
//  double nu_CMB_stat[static_spec_arr_size];
//  double nu_CMB_disk_RF_stat[static_spec_arr_size];
//  double n_CMB_stat[static_spec_arr_size];

 //   alloc_photons(&(pt_base->I_nu_CMB_stat),pt_base->core.nu_seed_size);
 //   alloc_photons(&(pt_base->I_nu_CMB_disk_RF_stat),pt_base->core.nu_seed_size);
 //   alloc_photons(&(pt_base->nu_CMB_stat),pt_base->core.nu_seed_size);
 //   alloc_photons(&(pt_base->nu_CMB_disk_RF_stat),pt_base->core.nu_seed_size);
 //   alloc_photons(&(pt_base->n_CMB_stat),pt_base->core.nu_seed_size);





//EC Disk
//  double nu_EC_Disk[static_spec_arr_size];
//  double nu_EC_Disk_obs[static_spec_arr_size];
//  double nuF_nu_EC_Disk_obs[static_spec_arr_size];
//  double nuFnu_EC_Disk_grid[static_spec_arr_size];

    alloc_photons(&(pt_base->Disk.ec.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->Disk.ec.spec.nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->Disk.ec.spec.nu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->Disk.ec.spec.nuFnu_obs),pt_base->core.nu_IC_size);


//  double nuFnu_Disk_grid[static_spec_arr_size];
//  double I_nu_Disk[static_spec_arr_size];
//  double J_nu_Disk_disk_RF[static_spec_arr_size];
//  double I_nu_Disk_disk_RF[static_spec_arr_size];
//  double nu_Disk[static_spec_arr_size];
//  double nu_Disk_obs[static_spec_arr_size];
//  double nu_Disk_disk_RF[static_spec_arr_size];
//  double nuF_nu_Disk_obs[static_spec_arr_size];
//  double n_Disk[static_spec_arr_size];

    alloc_photons(&(pt_base->Disk.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->Disk.spec.I_nu),pt_base->core.nu_seed_size);
    //alloc_photons(&(pt_base->J_nu_Disk_disk_RF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.I_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.nu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.nuFnu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->Disk.spec.n_nu),pt_base->core.nu_seed_size);



//EC BLR
//  double nuFnu_EC_BLR_grid[static_spec_arr_size];
//  double nuF_nu_EC_BLR_obs[static_spec_arr_size];
//  double nu_EC_BLR[static_spec_arr_size];
//  double nu_EC_BLR_obs[static_spec_arr_size];


    alloc_photons(&(pt_base->BLR.ec.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->BLR.ec.spec.nuFnu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->BLR.ec.spec.nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->BLR.ec.spec.nu_obs),pt_base->core.nu_IC_size);


//  double I_nu_BLR[static_spec_arr_size];
//  double nu_BLR[static_spec_arr_size];
//  double J_nu_BLR_disk_RF[static_spec_arr_size];
//  double nu_BLR_disk_RF[static_spec_arr_size];
//  double n_BLR[static_spec_arr_size];


    alloc_photons(&(pt_base->BLR.spec.I_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->BLR.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->BLR.spec.I_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->BLR.spec.nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->BLR.spec.n_nu),pt_base->core.nu_seed_size);



//EC DT
//  double nuFnu_EC_DT_grid[static_spec_arr_size];
//  double nu_EC_DT[static_spec_arr_size];
//  double nu_EC_DT_obs[static_spec_arr_size];
//  double nuF_nu_EC_DT_obs[static_spec_arr_size];

    alloc_photons(&(pt_base->DT.ec.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->DT.ec.spec.nu),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->DT.ec.spec.nu_obs),pt_base->core.nu_IC_size);
    alloc_photons(&(pt_base->DT.ec.spec.nuFnu_obs),pt_base->core.nu_IC_size);


//  double nuFnu_DT_grid[static_spec_arr_size];
//  double I_nu_DT[static_spec_arr_size];
//  double J_nu_DT_disk_RF[static_spec_arr_size];
//  double nu_DT_obs[static_spec_arr_size];
//  double nu_DT[static_spec_arr_size];
//  double nu_DT_disk_RF[static_spec_arr_size];
//  double n_DT[static_spec_arr_size];
//  double nuF_nu_DT_obs[static_spec_arr_size];

    alloc_photons(&(pt_base->DT.spec.nuFnu_grid),pt_base->core.nu_grid_size);

    alloc_photons(&(pt_base->DT.spec.I_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.I_nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.nu_obs),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.nu_DRF),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.n_nu),pt_base->core.nu_seed_size);
    alloc_photons(&(pt_base->DT.spec.nuFnu_obs),pt_base->core.nu_seed_size);

}

void alloc_photons(double ** pt,int size){
        printf("pre %p\n",*pt);
        unsigned int i;
        if (*pt){


        //   printf("freeing\n");

            //for (i = 0; i < size; i++) {

            //    printf("%e\n",*pt[i]);

            //}
            //free(*pt);
        }

        *pt= malloc(size * sizeof (double));
        printf("post %p\n",*pt);

    }

*/