//==============================================================================
//  FUNZIONI CHE COSTRUISCONO LA DISTRIBUZIONE ELETTRONICA
//==============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

/**
 * \file distrib_elettr.c
 * \author Andrea Tramacere
 * \date 19-09-2004
 * \brief funzioni per la
 * distribuzione energetica
 * degli elettroni sia nel caso
 * stazionario che per ET
 *
 */

//==============================================================
/// Genera la griglia su gamma che viene usata per integrare in gamma il
// il Sync ed l'IC. I punti con indice pari sono equispaziati nel logar
// tmo fra di loro, quelli con indice dispari sono la media fra il prece
// dente pari ed il successivo pari
//==============================================================

void Genera_griglia_gamma_N_log(struct spettro *pt, double * griglia_gamma_N_log) {
	unsigned long i;
    double delta_log;
    double log_a, log_b;
    if (pt->verbose>1) {
        printf("Generete log gamma_grid for N \n");
        printf("size is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    
    log_a = log10(pt->gmin_griglia);
    log_b = log10(pt->gmax_griglia);
    delta_log = (log_b - log_a) / ((double) pt->gamma_grid_size - 1);
    //PUNTI CON INDICE PARI LOG
    for (i = 0; i < pt->gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] = pow(10, (log_a + delta_log * (double) (i)));
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->griglia_gamma_Ne_log[i]);
    }
    //PUNTI CON INDICE DISPARI LIN
    for (i = 1; i < pt->gamma_grid_size; i += 2) {
        griglia_gamma_N_log[i] =
                (griglia_gamma_N_log[i - 1] + griglia_gamma_N_log[i + 1])*0.5;
        //printf("i=%d griglia_gamma_Ne_log=%e\n",i,pt->griglia_gamma_Ne_log[i]);
    }
}


//========================================
// Genera la  N[i] per e-
//========================================

void Genera_Ne(struct spettro *pt) {
    char *name;
    unsigned long i;

    if (pt->verbose>1) {
        printf("Set array per Ne \n");
        printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);
    }
    
    //Generate the gamma-grid for Ne(gamma)
    pt->griglia_gamma_Ne_log = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log);

    // Associa TIPO_DISTR  all stringa DISTR
    SetDistr(pt);

    //Fill Ne
    pt->Ne = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Fill_N(pt, pt->griglia_gamma_Ne_log, pt->Ne);
    
    
    
    //stationary frame
    pt->griglia_gamma_Ne_log_stat = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    for (i = 0; i < pt->gamma_grid_size; i++) {
		pt->griglia_gamma_Ne_log_stat[i]=pt->griglia_gamma_Ne_log[i]*pt->beam_obj;
	}
    // Associa TIPO_DISTR  all stringa DISTR
   

    //stationary frame
    pt->Ne_stat = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    for (i = 0; i < pt->gamma_grid_size; i ++) {
		pt->Ne_stat[i]=pt->Ne[i]*pt->beam_obj*pt->beam_obj;
	}
    
    
    //This flag is set to 1 to know that
    //N(gamma) has been properly initialized and filled
    pt->Distr_e_done = 1;
    pt->N_0e = pt->N_0;
	
	
	


    name = "distr-e.dat";
    Scrivi_N_file(pt, name, pt->griglia_gamma_Ne_log, pt->Ne);


}




//========================================
// Genera la  N[i] per pp ed e- secondari
//========================================

void Genera_Np_Ne_pp(struct spettro *pt) {
    char *name;
    double (*pf_distr) (struct spettro *, double x);
    printf("********** protons ***********\n");
    printf("set array for Np\n");
    printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);

    //Generate the gamma-grid for Np(gamma)
    pt->griglia_gamma_Np_log = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Np_log);

    // Associa TIPO_DISTR  all stringa DISTR
    SetDistr(pt);

    // Fill Np
    pt->Np = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Fill_N(pt, pt->griglia_gamma_Np_log, pt->Np);
    //This flag si set to 1 to know that
    //N(gamma) has been properly initialized and filled
    pt->Distr_p_done = 1;
    pt->N_0p = pt->N_0;

    name = "distr-p.dat";
    Scrivi_N_file(pt, name, pt->griglia_gamma_Np_log, pt->Np);


    // Secondaries e- from pp
    printf("****** secondary leptons *****\n");
    printf("set array for secondary Ne\n");
    printf("elements number is pt->gamma_grid_size=%d\n", pt->gamma_grid_size);

    //Generate the gamma-grid for Ne(gamma)
    pt->griglia_gamma_Ne_log = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Genera_griglia_gamma_N_log(pt, pt->griglia_gamma_Ne_log);

    //Set N to e- from pp
    double tmp_distr;
    tmp_distr = pt->TIPO_DISTR;
    pt->TIPO_DISTR = 8;

    //Fill Ne
    pt->Ne = (double*) malloc(pt->gamma_grid_size * sizeof (double));
    Fill_N(pt, pt->griglia_gamma_Ne_log, pt->Ne);
    //This flag si set to 1 to know that
    //N(gamma) has been properly initialized and filled
    pt->Distr_e_done = 1;
    pt->N_0e = pt->N_0;
    pf_distr = &N_distr_integranda;

    pt-> N_e_pp = integrale_trap_log_struct(pf_distr,
            pt,
            pt->gmin,
            pt->gmax,
            10000);
    printf("N_e_pp =%e\n", pt->N_e_pp);
    name = "distr-e-from-pp.dat";
    Scrivi_N_file(pt, name, pt->griglia_gamma_Ne_log, pt->Ne);
    pt->TIPO_DISTR = tmp_distr;
}





//========================================
// Sricve il file con  N[i]
//========================================

void Scrivi_N_file(struct spettro *pt, char *name, double *g, double *N) {
    double mass;
    char f_distr[512];
    FILE *fp_distr;
    sprintf(f_distr, "%s%s-%s", pt->path, pt->STEM, name);
    unsigned long i;
    
    fp_distr = fopen(f_distr, "w");
    if (fp_distr == NULL) {
        printf("warning non riesco ad aprire %s\n", f_distr);
        exit(1);
    }
    distr_e_header(fp_distr);

    if (pt->TIPO_DISTR != 6) {
        if (pt->TIPO_DISTR == 8 && strcmp(pt->PARTICLE, "hadrons") == 0) {
            mass = MEC2_TeV;
        }
        if (pt->TIPO_DISTR != 8 && strcmp(pt->PARTICLE, "hadrons") == 0) {
            mass = MPC2_TeV;
        }
        if (strcmp(pt->PARTICLE, "leptons") == 0) {
            mass = MEC2_TeV;
        }


        for (i = 0; i < pt->gamma_grid_size; i++) {
            if (N[i] > 0) {
                fprintf(fp_distr, "%e\t%e\t%e\t%e\t%e\t%e\n",
                        log10(g[i]),
                        log10(N[i]),
                        g[i],
                        N[i],
                        g[i] * mass,
                        N[i] / mass);
            }
        }
    }
    fclose(fp_distr);
}
//=====================================================





//========================================
// Trova il gmax da N[i]>o
//========================================

double Find_gmax(struct spettro *pt, double *N, double *g) {
	unsigned long i;
    double gmax;
    gmax = g[0];

    for (i = 0; i < pt->gamma_grid_size; i++) {
        if (N[i] > 0) {
            gmax = g[i];
        }
    }
    return gmax;
}
//=====================================================




//========================================
// RIEMPIE IL VETTORE  N[i]
//========================================

void Fill_N(struct spettro *pt, double * griglia_gamma_N_log, double * N) {
	unsigned long i, K;
    double a, gamma_piu, gamma_meno;
    double Npiu, Nmeno;
    double *N_File;
    double *G_File;
    int count, count_max;
    double GI, NI, NNI;
    //integranda Disre e
    double (*pf_norm) (struct spettro *, double x);
    char N_file[512];
    FILE *fp;


    //Fille N(gamma) array
    if (pt->TIPO_DISTR != 6) {
        //printf("***** Fill array of Particle distribution *****\n");



        //Normalization
        pt->N_0 = 1.0;
        //if distr is e- from pp no normalization to compute
        if (pt->Norm_distr == 1 && pt->TIPO_DISTR != 8 &&   pt->Norm_distr_L_e_Sync<=0.0) {
            pf_norm = &N_distr_integranda;
            //printf("Evaluate Normalization constant N_0\n");
            pt->N_0 = integrale_trap_log_struct(pf_norm, pt, pt->gmin, pt->gmax, 10000);
            //printf("N_0=%e\n", pt->N_0);

            for (i = 0; i < pt->gamma_grid_size; i++) {

            	N[i] = N_distr(pt, griglia_gamma_N_log[i]);
            	//    printf("i=%d griglia_gamma_Ne_log=%e Ne=%e\n",i,pt->griglia_gamma_Ne_log[i],N_distr(pt,pt->griglia_gamma_Ne_log[i]));
            }
        }





        if ( pt->Norm_distr_L_e_Sync>0.0) {
			pt->N=1.0;
			pt->N_0 = 1.0;

			for (i = 0; i < pt->gamma_grid_size; i++) {

				N[i] = N_distr(pt, griglia_gamma_N_log[i]);
				//    printf("i=%d griglia_gamma_Ne_log=%e Ne=%e\n",i,pt->griglia_gamma_Ne_log[i],N_distr(pt,pt->griglia_gamma_Ne_log[i]));
			}
			pt->Distr_e_done = 1;

			//printf ("PS=%e\n",Power_Sync_Electron(pt));
			pt->N=pt->Norm_distr_L_e_Sync/Power_Sync_Electron(pt);

			for (i = 0; i < pt->gamma_grid_size; i++) {

				N[i] = N_distr(pt, griglia_gamma_N_log[i]);
				//    printf("i=%d griglia_gamma_Ne_log=%e Ne=%e\n",i,pt->griglia_gamma_Ne_log[i],N_distr(pt,pt->griglia_gamma_Ne_log[i]));
			}

			//printf ("PS=%e\n",Power_Sync_Electron(pt));
			pt->Distr_e_done = 0;


       }

    }


    //=========================================
    // Leggo N(E) da file !!!EXPERIMENTAL
    //=========================================

    if (pt->TIPO_DISTR == 6) {
        printf("Leggo N(gamma) da file con formato gamma N(gamma) \n");
        sprintf(N_file, "%sN.dat", pt->path);
        fp = fopen(N_file, "r");
        count = 0;
        if (fp == NULL) {
            printf("Unable to open file %s \n ", N_file);
            exit(0);
        }
        while (!feof(fp)) {
            fscanf(fp, "%lf %lf \n", &GI, &NI);
            //printf("%e %e \n",GI,NI);
            //printf("count=%d\n",count);
            count++;
        }
        count_max = count - 1;
        printf("Elementi Vettore=%d  N_file[0:N-1=%d]\n", count_max + 1, count_max);
        fclose(fp);

        N_File = (double*) malloc((count_max + 1) * sizeof (double));
        G_File = (double*) malloc((count_max + 1) * sizeof (double));


        fp = fopen(N_file, "r");

        count = 0;
        while (!feof(fp) && count <= count_max) {
            fscanf(fp, "%lf %lf ", &GI, &NI);
            G_File[count] = GI;
            N_File[count] = NI;
            //printf("count=%d %e %e\n",count,GI,NI);
            count++;
        }
        count = 0;
        fclose(fp);

        printf("count_max=%d\n", count_max);
        //pt->gmin_griglia=G_File[0]/2;
        //if (pt->gmin_griglia<1.0){
        //    pt->gmin_griglia=1.0;
        //}
        //pt->gmax_griglia=G_File[count_max]*1.1;
        //pt->gmin=G_File[0];
        pt->gmax = G_File[count_max];

        printf("Aggirnamento parametri\n");
        printf("gmin=%e\n", pt->gmin);
        printf("gmax=%e\n", pt->gmax);
        printf("gmin_griglia=%e\n", pt->gmin_griglia);
        printf("gmax_griglia=%e\n", pt->gmax_griglia);




        //NORMALIZZO
        pt->N_0 = 1.0;
        if (pt->Norm_distr == 1) {
            pf_norm = &N_distr_integranda;
            printf("Evaluate Normalization constant N_0\n");
            //trapezoidal rule
            pt->N_0 = 0.0;
            for (i = 0; i < count_max; i++) {
                pt->N_0 += (G_File[i + 1] + G_File[i]) * N_File[i];
            }
            pt->N_0 *= 0.5;
            printf("Cost Norm=%e\n", pt->N_0);
        }

        //for(i=0;i<pt->gamma_grid_size;i++){
        //  printf("i=%d Gamma=%e\n",i,pt->griglia_gamma_Ne_log[i]);
        // }

        i = 0;
        N[0] = N_File[0];
        printf("i=%d Gamma=%e N_gamma=%e\n", i, griglia_gamma_N_log[i], N[i]);
        count = 0;

        //Interpolate The file distribution over the griglia_gamma grid
        for (i = 1; i < pt->gamma_grid_size; i++) {
            //printf("G=%e G_file=%e\n",pt->griglia_gamma_Ne_log[i],G_File[count]);
            if (griglia_gamma_N_log[i] < G_File[0] || griglia_gamma_N_log[i] > G_File[count_max]) {
                N[i] = 0;
            } else {
                //printf("<<i=%d G=%e G_file=%e\n",i,pt->griglia_gamma_Ne_log[i],G_File[count]);
                while (griglia_gamma_N_log[i] > G_File[count] && count < count_max - 1) {
                    //pt->griglia_gamma_Ne_log[i]=0;
                    //printf("   count=%d -G=%e G_file=%e\n",count,pt->griglia_gamma_Ne_log[i],G_File[count]);

                    count++;
                    //printf("   count=%d +G=%e G_file=%e\n",count,pt->griglia_gamma_Ne_log[i],G_File[count]);
                }
                //count--;
                //printf(">>i=%d count=%d  G=%e G_file=%e \n",i,count,pt->griglia_gamma_Ne_log[i],G_File[count]);

                //
                gamma_piu = G_File[count];
                gamma_meno = G_File[count - 1];
                //printf(">>interpolo fra %e %e \n",gamma_piu,gamma_meno);
                //
                Npiu = N_File[count];
                Nmeno = N_File[count - 1];

                a = ((griglia_gamma_N_log[i] - gamma_meno) / (gamma_piu - gamma_meno))*(Npiu - Nmeno);
                a += Nmeno;
                N[i] = a * pt->N / pt->N_0;
            }
            //printf("i=%d Gamma=%e N_gamma=%e\n\n",i,pt->griglia_gamma_Ne_log[i],pt->Ne[i]);
        }

        //LEGGI FILE FORMATO "GAMMA, N(GAMMA)differenzile"
        //CONVERI IN DOUBLE E NORMALIZZA
        //DETERMINA N_PUNTI G_MIN G_MAX
        //FAI IL CECK CON G_MIN G_MAX ANCHE PER GRIGLIA
        //EVENTUALMENTE CORREGGI E MANDA MESSAGGIO A VIEDO
        //QUESTO E' UN PUNTO CRITICO
        //ORA FAI COME PER IL CASO ET (VEDI SOPRA)
    }
}



//==============================================================
// N_distr
//==============================================================

double N_distr(struct spettro *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * questa funzione restituisce la distribuzione energetica                          \n
     * richiesti dalla funzione chiamante. Le funzioni chiamanti sono quelle per        \n
     * il calcolo degli spettri di sincrotrone ed di IC/EC                              \n
     *
     */
    double a, a1;
    double log_a;


    //Distribuzioni energetiche degli elettroni nel caso statico


    //PL
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 0) {
        return pt_N->N * pow(Gamma, -(pt_N->p)) / pt_N->N_0;
    }


    //PL CON EXP CUTOFF
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 1) {
        return pt_N->N * pow(Gamma, -(pt_N->p)) * exp(-Gamma / pt_N->gamma_cut) / pt_N->N_0;
    }

    //PL BROCKEN
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 2) {
        if (Gamma >= pt_N->gmin && Gamma <= pt_N->gamma_break)return pt_N->N * pow(Gamma, -pt_N->p) / pt_N->N_0;
        if (Gamma > pt_N->gamma_break && Gamma <= pt_N->gmax) {
            return pt_N->N * pow(pt_N->gamma_break, -(pt_N->p - pt_N->p_1)) * pow(Gamma, -pt_N->p_1) / pt_N->N_0;
        }
    }

    //LOG PARABOLA
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 3) {
        return pt_N->N * pow((Gamma / (pt_N->gamma0_log_parab)),
                (-pt_N->s - pt_N->r * log10(Gamma / pt_N->gamma0_log_parab))) / pt_N->N_0;
    }

    //LOG PARABOLA CON PICCO
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 4) {
        return pt_N->N * pow(10, (-pt_N->r * pow(log10(Gamma / pt_N->gammap_log_parab), 2))) / pt_N->N_0;
    }

    //LOG PARABOLA CON PL
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 5) {

        //Gamma<Gamma_c PL
        if (Gamma < pt_N->gamma0_log_parab) {
            return pt_N->N * pow((Gamma / pt_N->gamma0_log_parab),
                    (-pt_N->s)) / pt_N->N_0;
        }

        //Gamma>Gamma_c LOGPAR
        if (Gamma >= pt_N->gamma0_log_parab && Gamma <= pt_N->gmax) {
            return pt_N->N * pow((Gamma / (pt_N->gamma0_log_parab)),
                    (-pt_N->s - pt_N->r * log10(Gamma / pt_N->gamma0_log_parab))) / pt_N->N_0;
        }
    }

    //Spit
    double emin, norm, x, x_min;
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 7) {
        emin = 7.0 * pt_N->spit_cut;
        a = 0.;
        x = Gamma / pt_N->spit_cut;
        x_min = emin / pt_N->spit_cut;
        if (Gamma < emin) {
            a = (x) * exp(-x);
        } else {
            norm = (x_min) * exp(-x_min) / (pow(x_min, -pt_N->spit_index));
            //a1=1;
            a = (x) * exp(-x) + norm * pow(x, -pt_N->spit_index) * exp(-Gamma / pt_N->spit_cut1);
            //a*=min(1, exp((Gamma-emin*5)/emin*0.1));
            a *= pt_N->spit_ratio;
        }

        //printf('%e\n', pt_N->N*(a+a1)/pt_N->N_0);
        return pt_N->N * (a) / pt_N->N_0;
    }
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 8) {
        pt_N->Gamma = Gamma;
        //printf("gamma=%e rate=%e\n", pt_N->Gamma,rate_electrons_pp(pt_N, Gamma));

        return vluce_cm * pt_N->NH_pp * MEC2_TeV * bn_to_cm2 * rate_electrons_pp(pt_N, Gamma);
    }
    return 0.0;
}



//==============================================================
//   funzione integranda per la distribuzione degli e-
//    per calcolare il coeff di norm
//==============================================================

double N_distr_integranda(struct spettro *pt_N, double Gamma) {
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n
     * funzione che resitituisce le integrande per il calcolo  \n
     * del coefficiente di normalizzazione delle distribuzioni \n
     * elettroniche statiche                                   \n
     */

    double log_a, a, a1;
    //PL
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 0) {
        return pow(Gamma, -(pt_N->p));
    }

    //PL CON EXP CUTOFF
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 1) {
        return pow(Gamma, -(pt_N->p)) * exp(-Gamma / pt_N->gamma_cut);
    }

    //PL BROCKEN
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 2) {
        if (Gamma >= pt_N->gmin && Gamma <= pt_N->gamma_break)return pow(Gamma, -pt_N->p);
        if (Gamma > pt_N->gamma_break && Gamma <= pt_N->gmax) {
            return pow(pt_N->gamma_break, -(pt_N->p - pt_N->p_1)) * pow(Gamma, -pt_N->p_1);
        }
    }

    //LOG PARABOLA
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 3) {
        return pow((Gamma / (pt_N->gamma0_log_parab)),
                (-pt_N->s - pt_N->r * log10(Gamma / pt_N->gamma0_log_parab)));

    }


    //LOG PARABOLA CON PICCO
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 4) {
        return pow(10, (-pt_N->r * pow(log10(Gamma / pt_N->gammap_log_parab), 2)));
    }

    //LOG PARABOLA CON PL
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 5) {
        if (Gamma < pt_N->gamma0_log_parab) {
            return pow((Gamma / pt_N->gamma0_log_parab), (-pt_N->s));
        }
        if (Gamma >= pt_N->gamma0_log_parab && Gamma <= pt_N->gmax) {
            return pow((Gamma / (pt_N->gamma0_log_parab)),
                    (-pt_N->s - pt_N->r * log10(Gamma / pt_N->gamma0_log_parab)));
        }
        //return log_a;
    }

    //Spit
    double emin, norm, x, x_min;
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 7) {
        emin = 7.0 * pt_N->spit_cut;
        a = 0.;
        x = Gamma / pt_N->spit_cut;
        x_min = emin / pt_N->spit_cut;
        if (Gamma < emin) {
            a = (x) * exp(-x);

        } else {
            norm = (x_min) * exp(-x_min) / (pow(x_min, -pt_N->spit_index));
            a = (x) * exp(-x) + norm * pow(x, -pt_N->spit_index) * exp(-Gamma / pt_N->spit_cut1);
            a *= pt_N->spit_ratio;
        }



        //a=Gamma*exp(-Gamma/pt_N->spit_cut);
        //a1=pt_N->spit_ratio*pow(Gamma,-pt_N->spit_index);
        //a1*=min(1, exp(-(Gamma-pt_N->spit_cut1)/pt_N->spit_cut2) );
        //a1*=min(1, exp(Gamma-pt_N->spit_cut*2) );
        //exit(0);
        //printf("gamma=%e a+a1=%e\n",Gamma,a+a1);
        //a1=0;
        return a;
    }


    //Secondaris e Distribution has not analytical expression
    //it is taken from the N array, throug log-lin interpolation
    if (Gamma >= pt_N->gmin && Gamma <= pt_N->gmax && pt_N->TIPO_DISTR == 8) {

        return N_distr_interp(pt_N, Gamma, pt_N->griglia_gamma_Ne_log, pt_N->Ne);
    }

    return 0.0;
}

double N_distr_interp(struct spettro *pt, double Gamma, double *griglia_gamma, double *N) {
	unsigned long i;
    double gamma_piu, gamma_meno, Npiu, Nmeno, g, a;
    i = 0;
    while (griglia_gamma[i] < Gamma && i < pt->gamma_grid_size) {
        i++;
    }

    //printf("G=%e G_file=%e\n",pt->griglia_gamma_Ne_log[i],G_File[count]);
    if (i > 0 && i < pt->gamma_grid_size && N[i] > 0 && N[i - 1] > 0) {
        gamma_piu = log10(griglia_gamma[i]);
        gamma_meno = log10(griglia_gamma[i - 1]);
        Npiu = log10(N[i]);
        Nmeno = log10(N[i - 1]);
        g = log10(Gamma);
        a = ((g - gamma_meno) / (gamma_piu - gamma_meno))*(Npiu - Nmeno);
        a += Nmeno;
        //printf("%d %e %e %e %e %e %e\n",i,gamma_piu,gamma_meno,N[i],Npiu,Nmeno,a);
        return pow(10, a);
    } else {
        return 0;
    }
}
