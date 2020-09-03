#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"
/**
  * \file cosmo.c
  * \author Andrea Tramacere
  * \date 14-10-2004 \n
  * \brief file per distanza luminosita'
  */



double distanza_z(double z){
/**
   * \author Andrea Tramacere
   * \date 19-09-2004 \n
   * funzione integranda per il calcolo 
   * della distanza cosmologica per modelli 
   * a curvatura nulla (Omega_k=0);
   * 
   */
  //  double distanza;
  double a;
  a=pow((    (1+z)*(1+z)*(Omega_matter*z+1)  -z*(z+2)*Omega_lambda) ,-0.5);
  //printf("distanza=%e \n",a);
  return a;
}


double dist_lum_cm(double z_cosm){
  double dist;
  double (*pf) ( double);
  
  pf = &distanza_z;
  
  dist = (vluce_cm * 1.0e-5 / H_0)*(1.0 + z_cosm) * integrale_simp(pf, 0,z_cosm, 10000);
  dist *= parsec * 1.0e6 * 1.0e2;
  return dist;
}



double Distanza_Lum_analyt(double z){
  /**
   * \author Andrea Tramacere
   * \date 21-11-2004 \n
   * funzione per calcolare
   * la distanza luminosità
   * analiticamente per 
   * universo piatto 
   * formula presa da
   * THE ASTROPHYSICAL JOURNAL SUPPLEMENT SERIES, 120:49»50, 1999 January 1999 
   * ANALYTICAL FIT TO THE LUMINOSITY DISTANCE FOR FLAT COSMOLOGIES WITH A COSMOLOGICAL CONSTANT
   * author UE-LI PEN
   */
  double dl;
  dl=vluce_km/(H_0)*(1+z)*(eta_Distanza_Lum_analyt(1,Omega_matter)
			   -eta_Distanza_Lum_analyt((1/(1+z)),Omega_matter));
  return dl;
}

double eta_Distanza_Lum_analyt(double a, double Omega){
  double s,eta;
  s=pow((1-Omega)/Omega,1.0/3.0);
  eta=2*sqrt(s*s*s+1)*
    pow(( (1/pow(a,4.0))+
	  (-0.1540*(s/pow(a,3.0)))+
	  (+0.4304*(s*s/(a*a)))+
	  (+0.19097*(s*s*s/a))+
	  (+0.066941*s*s*s*s)
	  ),-(1.0/8.0)
	);
  //printf("eta=%e\n",eta);
  return eta;
}
