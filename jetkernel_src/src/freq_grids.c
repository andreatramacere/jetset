#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
//#include "libmia.h"
#include "Blazar_SED.h"

//========================
//GRID FUNCTIONS
//========================

void build_log_grid(double nu_start, double nu_stop, unsigned long SIZE, double * nu_grid){
	unsigned long I,I_MAX;
	double k,log_nu_start;
	k = (log10(nu_stop) - log10(nu_start));
	log_nu_start = log10(nu_start);
	I_MAX=SIZE-1;
	for (I = 0; I <=I_MAX; I++) {
		nu_grid[I] = pow(10, log_nu_start + k * (double) I / (double) I_MAX);
	}

}

unsigned long x_to_grid_index(double * nu_grid, double nu, unsigned long SIZE){
	unsigned long I=0;
	while (!(nu >= nu_grid[I] && nu <=nu_grid[I + 1] ) && (I<SIZE)){
		I++;
	}

	if (I==SIZE){
		I=-1;
	}

	return I;
}





//=========================================================================================






