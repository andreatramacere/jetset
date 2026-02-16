//=======================================================================
//
//   FUNZIONI MATEMATICHE
//=======================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"
#include "nrutil.h"

/**
 * \file funz_math.c
 * \author Andrea Tramacere
 * \date 27-04-2004
 * \brief funzioni matematiche
 *
 */


//double st_gamma(double x)
//{
//  return sqrt(2.0*pi/x)*pow(x/M_E, x);
//}


double st_gamma(double z)
{
  const int a = 12;
  static double c_space[12];
  static double *c = NULL;
  int k;
  double accm;

  if ( c == NULL ) {
    double k1_factrl = 1.0; /* (k - 1)!*(-1)^k with 0!==1*/
    c = c_space;
    c[0] = sqrt(2.0*M_PI);
    for(k=1; k < a; k++) {
      c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
	  k1_factrl *= -k;
    }
  }
  accm = c[0];
  for(k=1; k < a; k++) {
    accm += c[k] / ( z + k );
  }
  accm *= exp(-(z+a)) * pow(z+a, z+0.5); /* Gamma(z+1) */
  return accm/z;
}



//=====================================================
// funzione che resituisce la K di bessel di 5/3 per la F(x)
// UTILIZZA ROUTINES NUM. RECIPES.
//=====================================================

double bessel_K_53(struct blob *pt, double x) {
    double ri, rip, rkp, ord, a;
    ord = 5.0 / 3.0;
    bessik(x, ord, &ri, &a, &rip, &rkp);
    return a;
}

double bessel_K_23(struct blob *pt, double x) {
    double ri, rip, rkp, ord, a;
    ord = 2.0 / 3.0;
    bessik(x, ord, &ri, &a, &rip, &rkp);
    return a;
}


double bessel_K_pitch_ave(struct blob *pt, double x) {
    //The Astrophysical Journal, 334:L5-L8,1988 November 1
    double ri, rip, rkp, ord, a, K_43,K_13;
    ord = 4.0 / 3.0;
    bessik(x, ord, &ri, &K_43, &rip, &rkp);
    ord = 1.0 / 3.0;
    bessik(x, ord, &ri, &K_13, &rip, &rkp);
    a=(K_43*K_13 - 0.6*x*(K_43*K_43 - K_13*K_13));

    return a;
}

//=========================================================================================



//===================================================
// funzione che costruisce la tabelle delle funz. di Bessel
//===================================================

void tabella_Bessel(struct blob *pt_TB) {
    //double a, t_Bessel;

    double (*pf_53) (struct blob *, double x);
    unsigned int i;
    FILE *fp ;

    char f_bessel_file[static_file_name_max_legth];
    double _x,_y,_in_ave_y, _in_ave_x;
    /*** kernel Sync per alfa fisso e kernel delta ***/
        if (pt_TB->core.verbose > 0)
    {
        printf("Evaluation of Bessel Tables\n");
    }
	//printf("start=%e\n", pt_TB->Sync.x_Bessel_min = (pt_TB->Sync.spec.nu_min) /
	//	((3.0 * pt_TB->Sync.nu_B * (pt_TB->emitters.gmax_griglia * pt_TB->emitters.gmax_griglia) / 2.0) * pt_TB->Sync.sin_psi));
	//printf("stop=%e\n", pt_TB->Sync.x_Bessel_max = (pt_TB->Sync.spec.nu_max) /
	//	((3.0 * pt_TB->Sync.nu_B * (pt_TB->emitters.gmin_griglia * pt_TB->emitters.gmin_griglia) / 2.0) * pt_TB->Sync.sin_psi));
    //}
    pt_TB->Sync.x_Bessel_min = 1E-17;
    pt_TB->Sync.x_Bessel_max = 7.2E2;

    pt_TB->Sync.x_ave_Bessel_min = 1E-16;
    pt_TB->Sync.x_ave_Bessel_max = 3.5E2;

    pt_TB->Sync.log_x_Bessel_min = log10(pt_TB->Sync.x_Bessel_min);
    pt_TB->Sync.log_x_Bessel_max = log10(pt_TB->Sync.x_Bessel_max);
    

    pt_TB->Sync.log_x_ave_Bessel_min = log10(pt_TB->Sync.x_ave_Bessel_min);
    pt_TB->Sync.log_x_ave_Bessel_max = log10(pt_TB->Sync.x_ave_Bessel_max);

    //if(pt_TB->Sync.x_Bessel_max>Bessel_MAX)pt_TB->Sync.x_Bessel_max=Bessel_MAX;
    pf_53 = &bessel_K_53;
    //printf("SYSPATH =%s\n", pt_TB->core.SYSPATH);
    //return;
    //strcpy(f_bessel_file, pt_TB->core.SYSPATH);
    //printf('SYSPATH: %s\n', pt_TB->core.SYSPATH);
    //return;

    sprintf(f_bessel_file, "%s/F_Sync.dat",  pt_TB->core.SYSPATH);
    if (pt_TB->core.verbose>1) {
	//printf("gmax_griglia=%e -> x/xc_min=%e    gmin_griglia=%e -> x/xc_max=%e\n",
    //        pt_TB->emitters.gmax_griglia, pt_TB->Sync.x_Bessel_min,
     //       pt_TB->emitters.gmin_griglia, pt_TB->Sync.x_Bessel_max);
    	printf("Bessel Tables  in  file: %s\n", f_bessel_file);
    }


    fp = fopen(f_bessel_file, "r");

	build_log_grid( pt_TB->Sync.x_Bessel_min,  pt_TB->Sync.x_Bessel_max, static_bess_table_size,  pt_TB->Sync.F_Sync_x);
	build_log_grid( pt_TB->Sync.x_ave_Bessel_min,  pt_TB->Sync.x_ave_Bessel_max, static_bess_table_size, pt_TB->Sync.F_ave_Sync_x);

 
    //fclose(fp);

    fp = fopen(f_bessel_file, "w");
    printf("Bessel Tables  not found, was expected to be: %s\n", f_bessel_file);
    printf("Now GENERATING F_Sync.dat file\n");

    fprintf(fp, "# F_x F_y F_x_ave F_y_ave G_x G_y\n");
    
    for (i = 0; i < static_bess_table_size; i++) {
        //x = pt_TB->Sync.x_Bessel_min *
        //        pow((pt_TB->Sync.x_Bessel_max / pt_TB->Sync.x_Bessel_min),
        //        ((double) i) / ((double) elementi_tabelle - 1.0));
        
        //pt_TB->Sync.F_Sync_x[i] = x;
        pt_TB->Sync.F_Sync_y[i]= pt_TB->Sync.F_Sync_x[i] * integrale_trap_log_struct(pf_53, pt_TB, pt_TB->Sync.F_Sync_x[i], 1000, 1000);
        
        pt_TB->Sync.G_Sync_x[i]= pt_TB->Sync.F_Sync_x[i];
        pt_TB->Sync.G_Sync_y[i]= pt_TB->Sync.G_Sync_x[i] * bessel_K_23(pt_TB,pt_TB->Sync.G_Sync_x[i]);
        pt_TB->Sync.log_F_Sync_x[i] = log10(pt_TB->Sync.F_Sync_x[i]);
        pt_TB->Sync.log_G_Sync_x[i] = log10(pt_TB->Sync.G_Sync_x[i]);
        if (pt_TB->Sync.F_Sync_y[i]>0.0){
            pt_TB->Sync.log_F_Sync_y[i] = log10(pt_TB->Sync.F_Sync_y[i]);
        }
        else{
            pt_TB->Sync.log_F_Sync_y[i] = -300.0;
        }

        if ( pt_TB->Sync.G_Sync_y[i]>0.0){
            pt_TB->Sync.log_G_Sync_y[i] = log10(pt_TB->Sync.G_Sync_y[i]);
        }
        else{
            pt_TB->Sync.log_G_Sync_y[i] = -300.0;
        }

        //pt_TB->Sync.F_ave_Sync_x[i] = x;
        pt_TB->Sync.F_ave_Sync_y[i]= pt_TB->Sync.F_ave_Sync_x[i] *pt_TB->Sync.F_ave_Sync_x[i]*  bessel_K_pitch_ave(pt_TB,  pt_TB->Sync.F_ave_Sync_x[i]);
        pt_TB->Sync.log_F_ave_Sync_x[i] = log10(pt_TB->Sync.F_ave_Sync_x[i]);
        if(pt_TB->Sync.F_ave_Sync_y[i]>0.0){
            pt_TB->Sync.log_F_ave_Sync_y[i] = log10(pt_TB->Sync.F_ave_Sync_y[i]);
        }
        else{
            pt_TB->Sync.log_F_ave_Sync_y[i]=-300.0;
        }
        fprintf(fp, "%e %e %e %e %e %e\n",
                pt_TB->Sync.F_Sync_x[i] , pt_TB->Sync.F_Sync_y[i] , pt_TB->Sync.F_ave_Sync_x[i], pt_TB->Sync.F_ave_Sync_y[i], pt_TB->Sync.F_Sync_x[i], pt_TB->Sync.G_Sync_y[i]);
        
        //printf("i=%d i_max=%d x=%e F(x)=%e\n",i,elementi_tabelle,x,pt_TB->tabella_F[i][1]);
            
    }
    
    if (pt_TB->core.verbose > 0)
    {
        printf("i_max=%d elementi_tabelle=%d \n", i, static_bess_table_size);
        printf("F_Sync_x min=%e, %e\n", pt_TB->Sync.F_Sync_x[0], pt_TB->Sync.x_Bessel_min);
        printf("F_Sync_x max=%e, %e\n", pt_TB->Sync.F_Sync_x[static_bess_table_size - 1], pt_TB->Sync.x_Bessel_max);
        printf("F_Sync_y max=%e\n", pt_TB->Sync.F_Sync_y[static_bess_table_size - 1]);
    }
    fclose(fp);
    pt_TB->core.BESSEL_TABLE_DONE=1;
}
//=========================================================================================

//=====================================================================
//LOG-LINEAR INTERPOLATION
//=====================================================================


double log_lin_interp(double x,  double * x_grid, double x_min, double x_max, double *  y_grid , unsigned int SIZE, double emiss_lim){
	unsigned int ID;
	double y1,y2,x1,x2,a_c;
	ID=x_to_grid_index(x_grid,  x,   SIZE);

	if (ID<0 || ID>SIZE-2){
		return emiss_lim;
	}

	else if (x<x_min || x>x_max){
		return emiss_lim;
	}
    else if(y_grid[ID]<= emiss_lim || y_grid[ID+1]<= emiss_lim)
    {
        return emiss_lim;
    }
    
	else{
		y1 = log10(y_grid[ID]);
		y2 = log10(y_grid[ID + 1]);
		x1 = log10(x_grid[ID]);
		x2 = log10(x_grid[ID + 1]);
		a_c = (log10(x) - x1)*(y2 - y1) / (x2 - x1);
		a_c += y1;
		//printf("nu=%e, ID=%lu, SIZE=%d x1=%e x2=%e y1=%e y2=%e a_c=%e  return=%e    %e %e\n",nu,ID,SIZE,x1,x2,y1,y2,a_c,pow(10,a_c),flux_grid[ID],flux_grid[ID + 1]);

		return pow(10,a_c);
	}
}
//=========================================================================================

//=====================================================================
//LOG-LINEAR INTERPOLATION
//=====================================================================
//Lagrange Interpolation formula
//
double log_quad_interp(double x,  double * x_grid, double x_min, double x_max, double *  y_grid , unsigned int SIZE, double emiss_lim){
	unsigned int ID;
	double y1,y2,y3,x1,x2,x3;
	double a,b,c,f;
    ID=x_to_grid_index(x_grid,  x,   SIZE);

	if (ID<0 || ID>SIZE-2){
		return emiss_lim;
	}

	else if (x<x_min || x>x_max){
		return emiss_lim;
	}
    
    
	else{
        if (ID == SIZE-2){
            ID=ID-1;
        }
        if (y_grid[ID+2]<= emiss_lim && y_grid[ID+1]> emiss_lim){
            ID=ID-1;
        }
        //if (y_grid[ID+2]<= emiss_lim && y_grid[ID+1]<=  emiss_lim){
        //    ID=ID;
        //}
		y1 = log10(y_grid[ID]);
		y2 = log10(y_grid[ID + 1]);
        y3 = log10(y_grid[ID + 2]);
		x1 = log10(x_grid[ID]);
		x2 = log10(x_grid[ID + 1]);
        x3 = log10(x_grid[ID + 2]);
        x=log10(x); 
		a=(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3));
        b=(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3));
        c=(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2));
        f = y1*a +y2*b +y3*c;
        //printf("ID=%lu SIZE=%d x=%e x1=%e x2=%e x3=%e y1=%e y2=%e y3=%e \n",ID,SIZE,x,x1,x2,x3,y1,y2,y3);
		return pow(10,f);
	}
}
//=========================================================================================


//=====================================================================
//LOG-LOG INTERPOLATION
//=====================================================================
double log_log_interp(double log_x,  double * log_x_grid, double log_x_min, double log_x_max, double *  log_y_grid , unsigned int SIZE, double emiss_lim){
    unsigned int ID;
	double y1,y2,x1,x2,a_c;
	ID=x_to_grid_index(log_x_grid,  log_x,   SIZE);

	if (ID<0 || ID>SIZE-2){
		return emiss_lim;
	}
	else if (log_x<log_x_min || log_x>log_x_max){
			return emiss_lim;
		}
	else{
		y1 = log_y_grid[ID];
		y2 = log_y_grid[ID + 1];
		x1 = log_x_grid[ID];
		x2 = log_x_grid[ID + 1];
		a_c = (log_x - x1)*(y2 - y1) / (x2 - x1);
		a_c += y1;
		//printf("nu=%e, ID=%lu, SIZE=%d x1=%e x2=%e y1=%e y2=%e a_c=%e  return=%e    %e %e\n",log_x,ID,SIZE,x1,x2,y1,y2,a_c,pow(10,a_c),log_y_grid[ID],log_y_grid[ID + 1]);

		return pow(10,a_c);
	}
}
//=========================================================================================


//=====================================================================
//INTEGRAZIONE TRAPEZOIDALE CON INT APERTO E GRIGLIA  LINEARE CON ARRAYS
//=====================================================================

double trapzd_array_linear_grid(double *x, double *y, unsigned int SIZE)
{
    double I, y1, y2,delta_x;
    unsigned int INDEX;
    I = 0;
    //x1 = x[0];
    y1 = y[0];
    delta_x= x[1]-x[0];

    for (INDEX = 1; INDEX < SIZE; INDEX++)
    {
        y2 = y[INDEX];
        I += (y1 + y2);        
        y1 = y2;
        //printf("%e %e %d\n",nu2, Ptot, i);
    }
    return I * 0.5 * delta_x;
}

double trapzd_array_arbritary_grid( double *x, double *y, unsigned int SIZE)
{
    /**
     * \author Andrea Tramacere
     * \date 19-09-2004 \n

     */

    double I, y1, y2, x1, x2;
    unsigned int INDEX;

    I = 0;
    x1 = x[0];
    y1 =y[0];
    for (INDEX = 1; INDEX < SIZE; INDEX++)
    {
        x2 = x[INDEX];
        y2 = y[INDEX];
        I += (y1 + y2) * (x2 - x1);
        x1 = x2;
        y1 = y2;
        //printf("%e %e %d\n",nu2, Ptot, i);
    }
    return I * 0.5;
}
//=========================================================================================

//=========================================================
// INTEGRAZIONE ALLA SINPSON CON ARRAYS  E GRIGLIA EQUILOG CON MIDPOINT                
//=========================================================
double integr_simp_grid_equilog(double * x, double *y, unsigned int size) {
    double integr, delta;
    unsigned int ID;

    integr=0.;
    
    if (size % 2 == 0) {
         printf("grid size must be even");
         exit(0);
    }
    //this is necessary because you skip the mid point
    //in the loop when the grid is equilog spaced with midpoint
    for (ID = 1; ID < size - 1; ID=ID+2)
    {
       
        //QUESTO DELTA RIMANE QUI
        //PERCHE' LA GRIGLIA NON E' EQUISPACED
        //NON PUO ANDARE FUORI DAL LOOP
        delta=(x[ID+1]-x[ID-1]);
        integr+=(y[ID-1]+4.0*y[ID]+y[ID+1])*delta;
        //printf("%d ID=%d\n",ID);
       
    }
    //if (integr==0){
    //    printf("ID=%d size=%d \n",ID,size);
    //}
    return integr/6.0;

}

//============================================================================
// SIMPSON CON INT CHIUSO E GRIGLIA EQUI_LOG10 (log-spaced in gamma, integral in dγ)
//============================================================================


double integrale_simp_log_struct(double (*pf)(struct blob *, double),
                                 struct blob *pt,
                                 double a, double b,
                                 unsigned int n_intervalli)
{
    if (!(a > 0.0) || !(b > a)) return 0.0;

    // Simpson requires an even number of subintervals
    if (n_intervalli < 2) n_intervalli = 2;
    if (n_intervalli % 2 != 0) n_intervalli++;

    const double log10_a = log10(a);
    const double delta   = log10(b) - log10_a;
    const double denom   = (double)n_intervalli;

    // We will build points γ_i on a uniform log10 grid:
    // x_i = log10(a) + delta * (i / n_intervalli),  i=0..n_intervalli
    // γ_i = 10^{x_i}
    //
    // Composite Simpson with non-uniform steps in γ can be applied
    // per pair of intervals using the local widths:
    // For i=0,2,4,...:
    //   h0 = γ_{i+1} - γ_i
    //   h1 = γ_{i+2} - γ_{i+1}
    //   contribution = (h0+h1)/6 * [ f_i*(2 - h1/h0) + f_{i+1}*( (h0+h1)^2/(h0*h1) ) + f_{i+2}*(2 - h0/h1) ]
    //
    // This reduces to standard Simpson when h0==h1, and works well on log meshes.

    double integr = 0.0;

    // initialize i=0,1,2
    double x0 = log10_a;
    double x1 = log10_a + delta * (1.0 / denom);
    double x2 = log10_a + delta * (2.0 / denom);

    double g0 = a;
    double g1 = pow(10.0, x1);
    double g2 = pow(10.0, x2);

    double f0 = pf(pt, g0);
    double f1 = pf(pt, g1);
    double f2 = pf(pt, g2);

    for (unsigned int i = 0; i < n_intervalli; i += 2) {

        // enforce exact endpoint on the final point
        if (i + 2 == n_intervalli) {
            g2 = b;
            f2 = pf(pt, g2);
        }

        double h0 = g1 - g0;
        double h1 = g2 - g1;

        if (h0 > 0.0 && h1 > 0.0) {
            double r0 = h1 / h0;
            double r1 = h0 / h1;
            double H  = h0 + h1;

            // Composite Simpson for unequal adjacent steps (quadratic through 3 points)
            double term0 = f0 * (2.0 - r0);
            double term1 = f1 * (H * H / (h0 * h1));
            double term2 = f2 * (2.0 - r1);

            integr += (H / 6.0) * (term0 + term1 + term2);
        }

        // advance by 2: (0,1,2) -> (2,3,4)
        if (i + 2 >= n_intervalli) break;

        // new indices: i0=i+2, i1=i+3, i2=i+4
        g0 = g2; f0 = f2;

        double x3 = log10_a + delta * ((double)(i + 3) / denom);
        double x4 = log10_a + delta * ((double)(i + 4) / denom);

        g1 = pow(10.0, x3);
        g2 = pow(10.0, x4);

        f1 = pf(pt, g1);
        f2 = pf(pt, g2);
    }

    return integr;
}

//============================================================================
// INTEGRAZIONE TRAPEZOIDALE CON INT CHIUSO E GRIGLIA  EQUI_LOG
//============================================================================

double integrale_trap_log_struct(double (*pf) (struct blob *, double x), struct blob * pt, double a, double b, unsigned int n_intervalli) {
    double integr, k, griglia, ordinata, ordinata1, log10_a;
    double delta;
    integr = 0.0;
    griglia = 0.0;
    ordinata = 0.0;
    delta = log10(b) - log10(a);

    //CALCOLO in a e b fuori dall loop
    //per evitare errori di arrotondamento
    //sulla griglia logaritmica

    ordinata = a;
    log10_a=log10(a);
    ordinata1 = pow(10,log10_a+(delta * ((1) / ((double) n_intervalli - 1))));
    griglia = ordinata1 - ordinata;
    //printf("n_int=%e oridinata=%e ordinata1=%e,griglia=%e \n",n_intervalli,ordinata,ordinata1,griglia);
    integr += griglia * (pf(pt, ordinata) + pf(pt, ordinata1));

    ordinata = pow(10, log10_a+(delta * (((double) n_intervalli - 2) / ((double) n_intervalli - 1))));
    ordinata1 = b;
    griglia = ordinata1 - ordinata;
    //printf("n_int=%e oridinata=%e ordinata1=%e griglia=%e \n",n_intervalli,ordinata,ordinata1,griglia);
    integr += griglia * (pf(pt, ordinata) + pf(pt, ordinata1));


    for (k = 1; k < n_intervalli - 2; k = k + 1.0) {
        ordinata = pow(10, log10_a+(delta * (k / ((double) n_intervalli - 1))));
        ordinata1 = pow(10, log10_a+(delta * ((k + 1) / ((double) n_intervalli - 1))));
        griglia = ordinata1 - ordinata;
        integr += griglia * (pf(pt, ordinata) + pf(pt, ordinata1));
        //if(k==0 || k==n_intervalli-1) printf("n_int=%e k=%e,oridinata=%e,griglia=%e ba=%e\n",n_intervalli,k,ordinata,griglia,ba);
        //printf("integr=%e\n",integr);
    }
    //printf("integr=%e\n",integr);
    return 0.5 * integr;
}
//=========================================================================================




//=========================================================
// INTEGRAZIONE ALLA SINPSON CON INT CHIUSO E GRIGLIA LIN                
//=========================================================

double integrale_simp(double (*pf) ( double x), double a, double b, unsigned int n_intervalli) {
    double integr1, integr2, h, k;
    
    //CHECK on n_intervalli even
    if (n_intervalli % 2 != 0) {
        n_intervalli ++;
    }

    integr1 = 0.0;
    integr2 = 0.0;
    h = (b - a) / ((double) n_intervalli);
    for (k = 1.0; k < n_intervalli; k = k + 2.0) {
        integr1 += pf((a + (h * k)));
    }
    integr1 *= 4.0;
    //printf("intgr=%e\n",integr);
    for (k = 2.0; k < n_intervalli; k = k + 2.0) {
        integr2 += pf((a + (h * k)));
    }
    integr2 *= 2.0;
    //printf("intgr=%e x=%e \n",integr,a+(h*k));
    return h * (integr1 + integr2 + ((pf((a)) + pf((b))))) / 3.0;
}
//=========================================================================================





//=========================================================
// INTEGRAZIONE ALLA SINPSON CON INT CHIUSO E GRIGLIA LIN                
//=========================================================

double integrale_simp_struct(double (*pf) (struct blob *, double x), struct blob * pt, double a, double b, unsigned int n_intervalli) {
    double integr1, integr2, h, k;
    
    //CHECK on n_intervalli even
    if (n_intervalli % 2 != 0) {
        n_intervalli ++;
    }

    integr1 = 0.0;
    integr2 = 0.0;
    h = (b - a) / ((double) n_intervalli);
    for (k = 1.0; k < n_intervalli; k = k + 2.0) {
        integr1 += pf(pt, (a + (h * k)));
    }
    integr1 *= 4.0;
    //printf("intgr=%e\n",integr);
    for (k = 2.0; k < n_intervalli; k = k + 2.0) {
        integr2 += pf(pt, (a + (h * k)));
    }
    integr2 *= 2.0;
    //printf("intgr=%e x=%e \n",integr,a+(h*k));
    return h * (integr1 + integr2 + ((pf(pt, (a)) + pf(pt, (b))))) / 3.0;
}
//=========================================================================================

double theta_heaviside(double x){
    if (x>0.){
        return 1.0;
    }   else{
        return 0.;
    }
}

double V_region(struct blob *pt ) {
    double V;
    if (strcmp(pt->core.GEOMETRY, "spherical") == 0) {
	    V = four_by_three_pi * pt->core.R * pt->core.R* pt->core.R;
    }

	else if (strcmp(pt->core.GEOMETRY, "spherical_shell") == 0) {
        V = four_by_three_pi * pt->core.R_sh * pt->core.R_sh* pt->core.R_sh;
		V = four_by_three_pi * pt->core.R_ext_sh * pt->core.R_ext_sh * pt->core.R_ext_sh -V;
	}

	else {
		printf("GEOMETRY variable set to wrong value, posible spherical or spherical_shell \n");
		exit(0);
	}

    return V;
}

double S_sphere(struct blob *pt ) {
    double S;
    if (strcmp(pt->core.GEOMETRY, "spherical") == 0) {
	    S = four_pi * pt->core.R * pt->core.R;
    }

	else if (strcmp(pt->core.GEOMETRY, "spherical_shell") == 0) {
        S = four_pi * pt->core.R_ext_sh * pt->core.R_ext_sh;
	}

	else {
		printf("GEOMETRY variable set to wrong value, posible spherical or spherical_shell \n");
		exit(0);
	}

    return S;
}

double eval_beta_gamma(double gamma) {
    return sqrt(1 - 1 / (gamma * gamma));
}

double Larmor_radius(double gamma,double B, double sin_alpha){

	return (sqrt(gamma*gamma-1)*MEC2*sin_alpha)/(q_esu*B);
}

double Larmor_radius_to_gamma(double Larmor_radius,double B, double sin_alpha){
	double c=(q_esu*B*Larmor_radius)/(MEC2*sin_alpha);
	return sqrt(1+ c*c);
}

double get_beaming(double BulkFactor, double theta) {
	return 1.0 / (BulkFactor * (1 - eval_beta_gamma(BulkFactor) * cos(theta * Deg_to_Rad)));

}

//==============================================================
//FUNZIONE PER IL CALCOLO DELLA DERIVATA
// f'~(f(x+h)-f(x-h))/2
// h ~(machine_precision)^1/3*x
//==============================================================

double derivata(double (*pf) (struct blob *, double x), struct blob *pt_d, double x) {
    double h;
    h = x * 1e-7;
    if ((x - h) < pt_d->emitters.gmin) return (pf(pt_d, x + h) - pf(pt_d, x)) / (h);
    return (pf(pt_d, x + h) - pf(pt_d, x - h)) / (2 * h);
}
//=========================================================================================



//==========================================================
// FUNZIONE DI BESSEL MOD DA NUME RECIPES
//
//==========================================================

void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp) {

    double EPS;
    double FPMIN;
    int MAXIT;
    double XMIN;

    int i, l, nl;
    double a, a1,
            b, c,
            d, del, del1, delh, dels, e,
            f, fact, fact2, ff, gam1,
            gam2, gammi, gampl, h, p,
            pimu,
            q0,
            q1, q2, qnew, ril, ril1, rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xi, xi2, xmu, xmu2;

    EPS = 1.0e-10;
    FPMIN = 1.0e-30;
    MAXIT = 10000;
    XMIN = 2.0;

    if (x <= 0.0 || xnu < 0.0) {
    	printf("NRERROR bad arguments in bessik");
    	exit(0);
    }
    nl = (int) (xnu + 0.5);
    xmu = xnu - nl;
    xmu2 = xmu*xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    h = xnu*xi;
    if (h < FPMIN) h = FPMIN;
    b = xi2*xnu;
    d = 0.0;
    c = h;
    for (i = 1; i <= MAXIT; i++) {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;
        del = c*d;
        h = del*h;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (i > MAXIT){

    	printf("x too large in bessik; try asymptotic expansion");
    	exit(0);
    }
    ril = FPMIN;
    ripl = h*ril;
    ril1 = ril;
    rip1 = ripl;
    fact = xnu*xi;
    for (l = nl; l >= 1; l--) {
        ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    f = ripl / ril;
    if (x < XMIN) {
        x2 = 0.5 * x;
        pimu = pi*xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
        d = -log(x2);
        e = xmu*d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e) / e);
        beschb(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = exp(e);
        p = 0.5 * e / gampl;
        q0 = 0.5 / (e * gammi);
        c = 1.0;
        d = x2*x2;
        sum1 = p;
        for (i = 1; i <= MAXIT; i++) {
            ff = (i * ff + p + q0) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q0 /= (i + xmu);
            del = c*ff;
            sum += del;
            del1 = c * (p - i * ff);
            sum1 += del1;
            if (fabs(del) < fabs(sum) * EPS) break;
        }
        if (i > MAXIT){
        	printf("bessk series failed to converge");
        	exit(0);
        }
        rkmu = sum;
        rk1 = sum1*xi2;
    } else {
        b = 2.0 * (1.0 + x);
        d = 1.0 / b;
        h = delh = d;
        q1 = 0.0;
        q2 = 1.0;
        a1 = 0.25 - xmu2;
        q0 = c = a1;
        a = -a1;
        s = 1.0 + q0*delh;
        for (i = 2; i <= MAXIT; i++) {
            a -= 2 * (i - 1);
            c = -a * c / i;
            qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;
            q0 += c*qnew;
            b += 2.0;
            d = 1.0 / (b + a * d);
            delh = (b * d - 1.0) * delh;
            h += delh;
            dels = q0*delh;
            s += dels;
            if (fabs(dels / s) < EPS) break;
        }
        if (i > MAXIT){

        	printf("bessik: failure to converge in cf2");
        	exit(0);
        }
        h = a1*h;
        rkmu = sqrt(pi / (2.0 * x)) * exp(-x) / s;
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi;
    }
    rkmup = xmu * xi * rkmu - rk1;
    rimu = xi / (f * rkmu - rkmup);
    *ri = (rimu * ril1) / ril;
    *rip = (rimu * rip1) / ril;
    for (i = 1; i <= nl; i++) {
        rktemp = (xmu + i) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    *rk = rkmu;
    *rkp = xnu * xi * rkmu - rk1;
}

/* (C) Copr. 1986-92 Numerical Recipes Software !'K4$<%#110L&")|oV'4. */
//=========================================================================================




//==========================================================
//   FUNZIONI PER TEST NUMERIC
//
//==========================================================

double test_lunghezza_vettore(double mesh, double a, double b, int Max_elem) {
    /* se la mesh da un vettore troppo lungo */
    /* allora viene ricalcolata              */

    double Delta_log;
    int Num_int;
    Delta_log = log10(b) - log10(a);
    Num_int = (int) (Delta_log / mesh);
    printf("Num_int=%d\n", Num_int);
    if (Num_int > Max_elem) {
        printf("!!!!!!!!!!!!!!!!!!!!!Attenzione ho cambiato il parametro mesh da %e", mesh);
        mesh = (Delta_log / Max_elem);
        printf(" a%e\n", mesh);
    }
    return mesh;
}

double test_int(struct blob *pt_d, double x) {
    //  printf("f=%e, x=%e\n",x*x,x);
    return x*x;
}

double test_int1(double x) {
    // printf("f=%e, x=%e\n",x*x,x);
    return x*x;
}

double test_solid_anlge(struct blob *pt,double theta)
{
	//double sin_theta;
	//sin_theta = sqrt(1.0 - mu * mu);
	return 2 * pi;
}



//=========================================================================================


