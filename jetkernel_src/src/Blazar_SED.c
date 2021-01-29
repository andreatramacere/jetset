


// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <string.h>
// #include <unistd.h>
// //#include "libmia.h"
// #include "Blazar_SED.h"
// #define grad(x,t) (x>t) ? 0:1;
// #define ngrad(x,t) (x>t) ? 1:0;

// /**
// * \file Blazar_SED.c
// * \author Andrea Tramacere
// * \date 27-04-2004
// * \brief lettura input file
// * chimata sotto funzioni
// *
// */

// int main(int argc, char **argv) {

//    char help[] = "-h";
//    char help1[] = "-help";
//    char help2[] = "--help";
//    char verbose[] = "-verbose";
//    char tastiera[] = "tastiera";
//    char file[] = "file";
//    char make_log[] = "make_log";
//    char EVOLVE[] = "Evolve";

//    struct spettro spettro_root;
//    struct temp_ev temporal_root;
//    struct spettro *pt_spec;
//    struct temp_ev *pt_temporal;

//    pt_spec = &spettro_root;
//    pt_temporal = &temporal_root;

//    SYSPATH = getenv("BLAZARSED");
//    if (SYSPATH == NULL) {
//        printf("!!!Warning Blazar SED env variable not set\n");
//        printf("!!!please provide a path\n");
//        exit(0);
//    }
//    printf("BLAZARSED=%s\n", SYSPATH);

//    //================================================
//    //                 man page
//    //================================================
//    if (argc < 2 ||
//            strcmp(argv[1], help) == 0 ||
//            strcmp(argv[1], help1) == 0 ||
//            strcmp(argv[1], help2) == 0
//            ) {
//        manpage();
//        return 0;
//    }

//    printf("******************************************************************************\n");
//    printf("*         Blazar_SED_NEW               Version 3.1                             *\n");
//    printf("******************************************************************************\n\n");

//    printf("command lines args %d:\n", argc);
//    int i;
//    for (i = 0; i < argc; i++) {
//        printf("%s ", argv[i]);
//    }
//    printf("\n");

//    if (strcmp(argv[1], tastiera) == 0) {
//        printf("opzione non piu' supportata \n");
//        return 0;
//    }



//    if (strcmp(argv[1], file) == 0) {
//        printf("******************************************************************************\n");
//        printf("Source Parameters from File \n");
//        printf("******************************************************************************\n");
//        FileInput(argc, argv, pt_spec, pt_temporal);
//        printf("******************************************************************************\n");
//        printf("\n");

//        if (strcmp(argv[2], make_log) != 0) {
//            printf("******************************************************************************\n");
//            printf("Paramters Resume \n");
//            printf("******************************************************************************\n");
//            show_blob(spettro_root);
//            show_temp_ev(temporal_root);
//            printf("******************************************************************************\n");
//            printf("\n");

//            printf("******************************************************************************\n");
//            if (argc == 5 && strcmp(argv[4], verbose) == 0) {
//                printf("verbose mode \n");
//                spettro_root.verbose = 1;

//            } else {
//                spettro_root.verbose = 0;
//            }


//            if (argc <=5) {
//                printf("******************************************************************************\n");
//            printf("Run SED computation \n");
//            printf("******************************************************************************\n");
//                Init(pt_spec);
//                Run_SED(pt_spec);
//            }







//            if (argc > 5) {
//                if (strcmp(argv[5], EVOLVE) == 0) {
//                    pt_spec->Norm_distr = 0;
//                    temp_evolution(pt_spec, pt_temporal);
//                } else {
//                    Run_SED(pt_spec);
//                }
//            }
//            printf("******************************************************************************\n");
//            printf("\n");
//        }
//    }
//    return 0;
// }
// //=========================================================================================







// //======================================
// // Input from File
// //======================================

// void FileInput(int argc, char **argv, struct spettro *pt_spec, struct temp_ev *pt_temporal) {
//    FILE *fp_log, *fp_beaming_factor, *fp_PID_FILE, *fp_temporal;
//    double dato;
//    //struct spettro spettro_root;

//    char TIPO_DISTR[80];
//    char MODE[80], EV[80];

//    char stringa[80];
//    char stringa1[80];
//    char PID_FILE[100];
//    char make_log[] = "make_log";

//    // SECTIONS
//    char PROCESSES_AND_OUT_FILES[] = "PROCESSES_AND_OUT_FILES";
//    char SPECTRAL_AND_PHYSICAL_PARAMETERS[] = "SPECTRAL_AND_PHYSICAL_PARAMETERS";
//    char NUMERICAL_PARAMETERS[] = "NUMERICAL_PARAMETERS";
//    char EXTERNAL_COMPTON[] = "EXTERNAL_COMPTON";
//    char EVOLVE[] = "Evolve";
//    int i;
//    int out_log;
//    double test, prova;
//    double (*pf) (struct spettro *, double);


//    if (strcmp(argv[2], make_log) != 0) {
//        sprintf(pt_spec->path, argv[2]);
//        sprintf(pt_spec->STEM, argv[3]);
//        printf("PATH=%s\n", pt_spec->path);
//        printf("STEM=%s\n", pt_spec->STEM);
//        fp_log = fopen(argv[4], "r");
//        if (fp_log == NULL) {
//            perror("non trovo il file.log :-((  \n");
//            exit(1);
//        }
//        out_log = 0;
//        if (argc > 6 && strcmp(argv[5], EVOLVE) == 0) {
//            if (argc == 6) {
//                printf("temporal evolution input file name not provided\n");
//                manpage();
//                exit(1);
//            }
//            fp_temporal = fopen(argv[6], "r");
//            if (fp_temporal == NULL) {
//                perror("temporal evolution input file not found  \n");
//                exit(1);
//            }
//        }
//    }

//    if (strcmp(argv[2], make_log) == 0) {
//        printf("make ev input file\n");
//        fp_temporal = fopen(argv[4], "w");
//        out_log = 1;
//        printf("test %s\n", argv[4]);
//        //fclose(fp_temporal);
//    }

//    if (argc > 5) {
//        strcpy(EV, argv[5]);
//    }
//    if ((strcmp(EV, EVOLVE) == 0) || (out_log == 1)) {
//        //---------------Lettura Parametri da FIle---------------------------
//        printf(" Reading Temp Ev parameters\n");



//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf %s\n", &dato, stringa);
//            pt_spec->attesa_Sync_cooling = (int) dato;
//            printf("attesa_Sync_cooling %d\n", pt_spec->attesa_Sync_cooling);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "0/1    attesa_Sync_cooling\n");
//        }
//        printf("test\n");

//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf %s\n", &dato, stringa);
//            pt_spec->attesa_compton_cooling = (int) dato;
//            printf("attesa_compton_cooling %d\n", pt_spec->attesa_compton_cooling);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "0/1    attesa_compton_cooling\n");
//        }


//        //if (out_log == 0) {
//        //    fscanf(fp_temporal, "%lf %s\n", &dato, stringa);
//        //    pt_spec->attesa_EC_cooling = (int) dato;
//        //    printf("attesa_EC_cooling %d\n", pt_spec->attesa_EC_cooling);
//       // }
//        //if (out_log == 1) {
//        //    fprintf(fp_temporal, "0/1    attesa_EC_cooling\n");
//        //}


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->NUM_SET = (int) dato;
//            printf("NUM_SET=%d\n", pt_temporal->NUM_SET);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "int   NUM_SET(num_of_outfile)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->T_SIZE = (int) dato;
//            printf("T_SIZE=%d\n", pt_temporal->T_SIZE);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "int   T_SIZE(T_grid_size)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf %s\n", &dato, stringa);
//            pt_temporal->L_inj = dato;
//            printf("L_inj %d\n", pt_temporal->L_inj = dato);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    L_inj(erg/s)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->duration = dato;
//            printf("duration(s)=%e\n", pt_temporal->duration);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    duration(R/C)\n");
//        }

//        //fscanf(fp,"%lf  %s",&dato,stringa);
//        //c2=dato;
//        //printf("T_ACC(R/C)=%e\n",c2);

//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->t_D0 = dato;
//            printf("t_D=%e (s)",pt_temporal->t_D0);
//            printf("Diff_Coeff=1/t_D=%e (1/1s) \n", 1.0 / pt_temporal->t_D0);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    t_D=1/Diff_Coeff(s)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->Diff_Index = dato;
//            printf("Diff_Index =%e\n", pt_temporal->Diff_Index);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    Diff_Index\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->t_A0 = dato;
//            printf("t_A=%e(/s) \n",  pt_temporal->t_A0);
//            printf("Acc_Coeff=1/t_A=%e(1/s) \n",  1.0/pt_temporal->t_A0);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    1/Acc_Coeff(s)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->Acc_Index = dato;
//            printf("Acc_Index =%e\n", pt_temporal->Acc_Index);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    Acc_Index\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->T_esc_Coeff = dato;
//            printf("T_ESC_COEFF(R/c)_DOUB =%e\n", pt_temporal->T_esc_Coeff);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    T_ESC_COEFF(R/c)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->Esc_Index = dato;
//            printf("ESC_Index =%e\n", pt_temporal->Esc_Index);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    ESC_index\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->TStart_Inj = dato;
//            printf("T_Start_Inj(s)_DOUB =%e\n", pt_temporal->TStart_Inj);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    T_Start_Inj(s)\n");
//        }



//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->TStop_Inj = dato;
//            printf("T_Stop_Inj(s)_DOUB =%e\n", pt_temporal->TStop_Inj);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    T_Stop_Inj(s)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->TStart_Acc = dato;
//            printf("T_Start_Acc(s)_DOUB =%e\n", pt_temporal->TStart_Acc);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    T_Start_Acc(s)\n");
//        }


//        if (out_log == 0) {
//            fscanf(fp_temporal, "%lf  %s", &dato, stringa);
//            pt_temporal->TStop_Acc = dato;
//            printf("T_Stop_Acc(s)_DOUB =%e\n", pt_temporal->TStop_Acc);

//            fclose(fp_temporal);
//        }
//        if (out_log == 1) {
//            fprintf(fp_temporal, "double    T_Stop_Acc(s)\n");
//            fclose(fp_temporal);
//        }


//    }





//    if (strcmp(argv[2], make_log) == 0) {
//        printf("make spec input file\n");
//        fp_log = fopen(argv[3], "w");
//        out_log = 1;
//    }

//    if (out_log == 0) {
//        printf("******************************************************************************\n");
//        printf("*         Lettura dei parametri dal file di inizializzazione                 *\n");
//        printf("******************************************************************************\n");
//    }

//    if (out_log == 1) {
//        printf("******************************************************************************\n");
//        printf("*         Scrivo il modello del file di inizializzazione                     *\n");
//        printf("******************************************************************************\n");
//    }

//    while (feof(fp_log) != 1 || (out_log == 1)) {

//        // SECTION PROCESSES_AND_OUT_FILES
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", stringa1, stringa);
//            printf("************ SECTION %s ************************\n", stringa);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "**************SECTION    %s\n", PROCESSES_AND_OUT_FILES);
//        }

//        // Particles
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", stringa, stringa1);
//            sprintf(pt_spec->PARTICLE, stringa);
//            printf("PARTICLES %s\n", pt_spec->PARTICLE);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "leptons \t\t\t PARTICLES(leptons,hadrons)\n");
//        }


//        // Synch
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->do_Sync = (int) dato;
//            printf("do_Sync %d\n", pt_spec->do_Sync);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1 \t\t\t do_Sync(0=no/1=yes/2=self_abs)\n");
//        }

//        // compton
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->do_SSC = (int) dato;
//            printf("do_SSC %d\n", pt_spec->do_SSC);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1 \t\t\t do_SSC(0=no/1=si)\n");
//        }

//        // mode
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", MODE, stringa);
//            sprintf(pt_spec->MODE, MODE);
//            printf("mode %s\n", MODE);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "fast \t\t\t mode(fast/accurare)\n");
//        }

//        // SECTION SPECTRAL_AND_PHYSICAL_PARAMETERS
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", stringa1, stringa);
//            printf("*****************SECTION %s*****************\n", stringa);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "**************SECTION    %s\n", SPECTRAL_AND_PHYSICAL_PARAMETERS);
//        }



//        // gamma min griglia
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gmin_griglia = dato;
//            printf("gmin_griglia %e\n", pt_spec->gmin_griglia);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1 \t\t\t gmin_griglia(double)\n");
//        }

//        // gamma min griglia
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gmax_griglia = dato;
//            printf("gmax_griglia %e\n", pt_spec->gmax_griglia);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1e5 \t\t\t gmax_griglia(double)\n");
//        }


//        // nu_seed_size
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_seed_size = (int) dato; /* determina quanti punti calcolare per il Sync */
//            printf("nu_seed_size %d\n", pt_spec->nu_seed_size);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "100 \t\t\t mesh_spettro_Sync(int)\n");
//        }

//        // mesh_spettro_comp
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_IC_size = (int) dato; /* determina quanti punti calcolare per il compton */
//            printf("nu_IC_size %d\n", pt_spec->nu_IC_size);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "200 \t\t\t mesh_spettro_comp(int)\n");
//        }

//        // nu_start_Sync
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_start_Sync = dato;
//            printf("nu_start_Sync %e Hz\n", pt_spec->nu_start_Sync);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e12 \t\t\t nu_start_Sync(double)\n");
//        }

//        // nu_stop_Sync
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_stop_Sync = dato;
//            printf("nu_stop_Sync %e Hz\n", pt_spec->nu_stop_Sync);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e20 \t\t\t nu_stop_Sync(double)\n");
//        }

//        // nu_start_comp
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_start_SSC = dato;
//            printf("nu_start_comp %e Hz\n", pt_spec->nu_start_SSC);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e18 \t\t\t nu_start_compton(double)\n");
//        }

//        // nu_stop_comp
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_stop_SSC = dato;
//            printf("nu_stop_comp %e Hz\n", pt_spec->nu_stop_SSC);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e27 \t\t\t nu_stop_compton(double)\n");
//        }


//        // grid somma mesh
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_sum_size = (int) dato; // determina quanti punti calcolare per la somma c+s
//            printf("mesh_somm %d\n", pt_spec->nu_sum_size);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "200 \t\t\t mesh_spettro_somma(int)\n");
//        }

//        //nu_start_grid
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_start_grid = dato;
//            printf("nu_start_grid %e Hz\n", pt_spec->nu_start_grid);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e9 \t\t\t nu_start_grid(double)\n");
//        }

//        //nu_stop_grid
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_stop_grid = dato;
//            printf("nu_stop_grid %e Hz\n", pt_spec->nu_stop_grid);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e20 \t\t\t nu_stop_grid(double)\n");
//        }

//        // B campo mag
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->B = dato;
//            printf("B %e G\n", pt_spec->B);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "0.1 \t\t\t B_campo_mag_in_Gauss(double)\n");
//        }

//        // R raggio blob
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->R = pow(10.0, dato); /* raggio blob o raggio shell cilindrica in centimetri */
//            printf("R %ecm\n", pt_spec->R);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "15 \t\t\t log10(R_blob(cm,double))\n");
//        }

//        // ALTRI DATI
//        // fattore di beaming
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->BulkFactor = dato;
//            printf("Bulk Factor=%e\n", pt_spec->BulkFactor);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "10 \t\t\t Gamma(BulkFactor,double)\n");
//        }


//        // viewing angle
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->theta = dato;
//            printf("viewing angle=%e\n", pt_spec->theta);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0 \t\t\t viewing_angle\n");
//        }



//        // red shifht
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->z_cosm = dato;
//            printf("red shift %e\n", pt_spec->z_cosm);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "0.1 \t\t\t z(red_shift,double)\n");
//        }


//        // NH_pp
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->NH_pp = dato;
//            printf("NH_pp %e\n", pt_spec->NH_pp);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0 \t\t\t  NH_pp(#cm^-3,double)\n");
//        }


//        // Norm Distr
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->Norm_distr = (int) dato; // 1=Nomr, 0 not Norm
//            printf("NormDistr %d  \n", pt_spec->Norm_distr);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1  \t\t\t 1=Distribution_Normalized,0=Not\n");
//        }



//        // N elettroni
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->N = dato; /*numero di elettroni per cm^3 */
//            printf("N %e e/cm^3\n", pt_spec->N);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "10 \t\t\t N_e-_#e-/cm^3(double)\n");
//        }

//        // Tipo di distribuzione Stazionaria
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", TIPO_DISTR, stringa);
//            sprintf(pt_spec->DISTR, TIPO_DISTR);
//            printf("Tipo di distribuzione %s\n", TIPO_DISTR);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "pl \t\t\t TIPO_DISTR_STAZ:pl;plc;bkn;lp;lpEp;lppl\n");
//        }

//        // p indice_pl
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->p = dato;
//            printf("p %e indic. distr_e- power law\n", pt_spec->p);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "2.0 \t\t\t p(indice_pl,double)\n");
//        }

//        // p_1 indice per broken pl
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->p_1 = dato;
//            printf("p_1 %e indic. distr_e- power law broken \n", pt_spec->p_1);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "3.0 \t\t\t p1(indice_pl_broken,double)\n");
//        }

//        // gamma_break
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gamma_break = pow(10.0, dato);
//            printf("gamma_break %e\n", pt_spec->gamma_break);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "4.0 \t\t\t log10(gamma_breack_distr_e-(double))\n");
//        }

//        // gamma cut off distr e- static
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gamma_cut = pow(10.0, dato);
//            printf("gamma_cut %e\n", pt_spec->gamma_cut);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "4.0 \t\t\t log10(gamma_cut_off_distr_e-(double))\n");
//        }


//        //Spitkovsky
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->spit_index = dato;
//            printf("spit_index %e\n", pt_spec->spit_index);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "2.4 \t\t\t spit_index\n");
//        }

//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->spit_ratio = dato;
//            printf("spit_ratio %e\n", pt_spec->spit_ratio);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "2.4 \t\t\t spit_ratio\n");
//        }

//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->spit_cut = pow(10.0, dato);
//            printf("spit_cut %e\n", pt_spec->spit_cut);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "4 \t\t\t lgo10(spit_cut)\n");
//        }

//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->spit_cut1 = pow(10.0, dato);
//            printf("spit_cut1 %e\n", pt_spec->spit_cut1);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "6 \t\t\t lgo10(spit_cut1)\n");
//        }

//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->spit_cut2 = pow(10.0, dato);
//            printf("spit_cut2 %e\n", pt_spec->spit_cut2);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "6.5 \t\t\t lgo10(spit_cut2)\n");
//        }




//        // indice r,s distribuzioni Log par

//        // r   indice r log par LP+LPPL
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->r = dato;
//            printf("r   indic.  r distr e- log par %e\n", pt_spec->r);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.5 \t\t\t r(indice_log_par,double)\n");
//        }

//        // s   indice s log par LPPL
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->s = dato;
//            printf("s   indic.  s distr e- log par  %e\n", pt_spec->s);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "2.0 \t\t\t s(indice_log_par,double)\n");
//        }

//        // indice Gamma_0 distribuzioni Log par
//        // gamma_0   indice gamma_0 log par statica LP+LPPL
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gamma0_log_parab = pow(10.0, dato);
//            printf("gamma0_log_parab   indic.  gamma0  distr e- log par LP,LPPL %e\n", pt_spec->gamma0_log_parab);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "4.0 \t\t\t log10(gamma0_log_par(indice_gamma0_log_pa-LP-LPPL,double))\n");
//        }

//        // indice Gammap distribuzioni Log par
//        // gamma_p   indice gamma_p log par statica LPEP
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gammap_log_parab = pow(10.0, dato);
//            printf("gammap_log_parab   indic.  gamma_p  distr e- log par LPEP %e\n", pt_spec->gammap_log_parab);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "4.0 \t\t\t log10(gammap_log_par(indice_gammap_log_par-LPEP,double))\n");
//        }

//        // Parametri Gamma griglia

//        // gamma min distr e- statico/inj
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gmin = pow(10.0, dato);
//            printf("gamma_min %e\n", pt_spec->gmin);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0 \t\t\t log10(gamma_min_distr_e-(double))\n");
//        }

//        // gamma max distr e- statico/inj*
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->gmax = pow(10.0, dato);
//            printf("gamma_max %e\n", pt_spec->gmax);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "6.0 \t\t\t log10(gamma_max_distr_e-(double))\n");
//        }

//        // SECTION EXTERNAL COMPTON
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", stringa1, stringa);
//            printf("*****************SECTION %s*****************\n", stringa);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "**************SECTION    %s\n", EXTERNAL_COMPTON);
//        }

//        // EC Disk
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->do_EC_Disk = (int) dato;
//            printf("do_EC_Disk %d\n", pt_spec->do_EC_Disk);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1 \t\t\t do_EC_Disk(0/1)\n");
//        }

//        // EC BLR
// 		if (out_log == 0) {
// 			fscanf(fp_log, "%lf %s\n", &dato, stringa);
// 			pt_spec->do_EC_BLR = (int) dato;
// 			printf("do_EC_BLR %d\n", pt_spec->do_EC_BLR);
// 		}
// 		if (out_log == 1) {
// 			fprintf(fp_log, "1 \t\t\t do_EC_BLR(0/1)\n");
// 		}
// 		 // EC DT
// 		if (out_log == 0) {
// 			fscanf(fp_log, "%lf %s\n", &dato, stringa);
// 			pt_spec->do_EC_DT = (int) dato;
// 			printf("do_EC_DT %d\n", pt_spec->do_EC_DT);
// 		}
// 		if (out_log == 1) {
// 			fprintf(fp_log, "1 \t\t\t do_EC_DT(0/1)\n");
// 		}

//        // Diks Type
//        if (out_log == 0) {
//            fscanf(fp_log, "%s %s\n", TIPO_DISTR, stringa);
//            sprintf(pt_spec->disk_type, TIPO_DISTR);
//            printf("Disk Type %s \n", pt_spec->disk_type);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "BB \t\t\t Disk_Type(BB,MultiBB,Mono)\n");
//        }

//        // nu_start_EC_BLR
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_start_EC_BLR = dato;
//            printf("nu_start_EC_BLR %e Hz\n", pt_spec->nu_start_EC_BLR);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e18 \t\t\t nu_start_BLR_EC(double)\n");
//        }

//        // nu_stop_BLR_EC
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_stop_EC_BLR = dato;
//            printf("nu_stop_EC_BLR %e Hz\n", pt_spec->nu_stop_EC_BLR);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e27 \t\t\t nu_stop_EC_BLR(double)\n");
//        }

//        // R_inner
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->R_inner_Sw = dato;
//            printf("R_inner_Sw %e(R_s)\n", pt_spec->R_inner_Sw);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "3.0 \t\t\t R_inner_Sw(Rs)\n");
//        }

//        // R_ext
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->R_ext_Sw = dato;
//            printf("R_ext_Sw %e (Rs)\n", pt_spec->R_ext_Sw);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "500.0 \t\t\t R_ext_Sw(Rs)\n");
//        }

//        // Lum Disk
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->L_Disk = pow(10.0, dato);
//            printf("Lum Diks %e\n", pt_spec->L_Disk);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "48 \t\t\t log10(L_Disk)\n");
//        }

//        // accr eff
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->accr_eff = dato;
//            printf("accr eff %e\n", pt_spec->accr_eff);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "0.1 \t\t\t log10(accr_eff)\n");
//        }


//        // tau BLR
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->tau_BLR = dato;
//            printf("tau BLR %e\n", pt_spec->tau_BLR);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1e-3 \t\t\t tau_BLR\n");
//        }

//        // T disk
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            //pt_spec->T_disk=dato;
//            pt_spec->T_Disk = dato;
//            printf("T disk max%e (K)\n", pt_spec->T_Disk);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1e5 \t\t\t T_Disk(K)(T_max_for_MultiBB)\n");
//        }

//        // Dist disk BLR
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->R_BLR_in = pow(10, dato);
//            printf("dist disk BLR %e\n", pt_spec->R_BLR_in);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "18 \t\t\t log10(R_BLR_in)(cm)\n");
//            //return spettro_root;
//        }

//        // nu_start_DT_EC
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_start_EC_DT = dato;
//            printf("nu_start_EC_DT %e Hz\n", pt_spec->nu_start_EC_DT);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e13 \t\t\t nu_start_EC_DT(double)\n");
//        }

//        // nu_stop_EC_DT
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->nu_stop_EC_DT = dato;
//            printf("nu_stop_EC_DT %e Hz\n", pt_spec->nu_stop_EC_DT);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1.0e27 \t\t\t nu_stop_EC_DT(double)\n");
//        }

//        // T Dusty Torus
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            //pt_spec->T_disk=dato;
//            pt_spec->T_DT = dato;
//            printf("T Dusty Torus %e\n", pt_spec->T_DT);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1e2 \t\t\t T_DT(K)(T_Dusty_Torus)\n");
//        }

//        // Dist disk DT
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->R_DT = pow(10, dato);
//            printf("dist disk DT %e\n", pt_spec->R_DT);
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "18.5 \t\t\t log10(R_DT)(cm)\n");
//            //return spettro_root;
//        }

//        // tau DT
//        if (out_log == 0) {
//            fscanf(fp_log, "%lf %s\n", &dato, stringa);
//            pt_spec->tau_DT = dato;
//            printf("tau BLR %e\n", pt_spec->tau_DT);
//            fclose(fp_log);
//            return;
//        }
//        if (out_log == 1) {
//            fprintf(fp_log, "1e-3 \t\t\t tau_DT\n");
//            fclose(fp_log);

//            return;
//        }


//    }
//    return;
//    printf("******************************************************************************\n");
//    printf("******************************************************************************\n\n");

// }



// //======================================
// // Man page
// //======================================

// void manpage() {
//    char help[] = "\n\
//    Blazar_SED man page:\n\n\
//    -to run:\n\
//    \tBlazar_SED_NEW file output_dir_path stem input_file \n\
//    \tBlazar_SED_NEW file output_dir_path stem spec_input_file Evolve ev_input_file\n\
//    \texample:\n\
//    \tBlazar_SED_NEW file /home/andrea/sim/sim1 001 input.log\n\n\
//    \ttempoaral evolution:\n\
//    \tBlazar_SED_NEW file /home/andrea/sim/sim1 001 spec_input.log Evolve ev_input.log\n\n\
//    -for help:\n\
//    \tBlazar_SED_NEW \n\
//    \tBlazar_SED_ASCD -h \n\
//    \tBlazar_SED_NEW --help\
//    \n\n\
//    -to generate a template input file:\n\
//    \tBlazar_SED_NEW file make_log spec_input_file_name ev_input_file_name\
//    \n\n\n";
//    printf("%s", help);
//    return;
// }

