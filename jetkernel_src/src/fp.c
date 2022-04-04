#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Blazar_SED.h"

double Inj_temp_prof(double t,struct temp_ev *pt){

	if (t>pt->TStart_Inj && t<pt->TStop_Inj) {
            //return 1.0+pt->Inj_temp_slope*((t-pt->TStart_Inj)/(pt->TStop_Inj-pt->TStart_Inj));
            //if (fabs(t-pt->TStart_Inj)<100 || fabs(t-pt->TStop_Inj)<100)
            //    {
                //printf('%e %e\n',fabs(t-pt->TStart_Inj),fabs(t-pt->TStop_Inj));
            //    return 1.0;
           // }

            //else return 0.0;
        return 1.0;
	}
	else return 0.0;
}

double f_Dp(double gamma,struct temp_ev *pt){
	if (gamma<pt->Gamma_Max_Turb_L_coher){
		return pt->Diff_Coeff*pow(gamma,pt->Diff_Index);
	}
	else if((gamma>=(pt->Gamma_Max_Turb_L_coher)) && gamma<((pt->Gamma_Max_Turb_L_max))){
		return pt->Diff_coeff_CD;
	}
	else{
		//return pt->Diff_coeff_CD;
		return 1E60;
	}


}

double f_DAp(double gamma,struct temp_ev *pt){
	if (gamma<pt->Gamma_Max_Turb_L_coher){
		return pt->Diff_Coeff*2*pow(gamma,pt->Diff_Index-1);
	}
	else if((gamma>=(pt->Gamma_Max_Turb_L_coher)) && gamma<((pt->Gamma_Max_Turb_L_max))){
		return pt->Diff_coeff_CA*pow(gamma,-1);
	}
	else{
		return 1E60;


	}
}

double f_A(double gamma,struct temp_ev *pt){
	return pt->Acc_Coeff*pow(gamma,pt->Acc_Index);
}



//---------------------------------------------------------------------
double f_Tesc(double x, double coeff, double index){
	//Definire il tempo di fuga
	//gamma=x+1
	return    coeff*pow(x,index);
	//return pt->T_esc_Coeff;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
double Cooling(double x,struct temp_ev *pt, struct blob *pt_spec){
	double cooling;
	//gamma-1=x
	//gamma=x+1
	// printf("Sync cooling =%e\n",Sync_cool(pt_spec->B,x+1));

	cooling=0;
	if ((x+1)>1){
		if (pt->do_Sync_cooling>0){
			//printf('ciccio s \n');
			cooling=Sync_cool(pt_spec,x+1);

		}

		if (pt->do_Compton_cooling>0){
			//printf('ciccio c  \n');
			cooling+=compton_cooling(pt_spec,pt,x+1);
		}
		if(pt->do_Expansion>0 && pt->t>pt->t_jet_exp && pt->do_Adiabatic_cooling>0){
			cooling+=(x+1)/(Adiabatic_Cooling_time(pt, pt_spec,pt->R_jet_t));
		}
	}
	return cooling;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
//Termine Diffusivo-> Stocastico
double Cfp(double x, struct temp_ev *pt){
	return f_Dp(x+1,pt);
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
//Accelerazione Sistematica
double f_Acc(double x, struct temp_ev *pt){
	double D_A,A;
	//D_A=2*D/x=2*D_0*X^(q-1)
	//A=A_0*x^(Acc_Index)
	//x=gamma-1
	D_A=f_DAp(x+1,pt);
	A=f_A(x+1,pt);

	return A+D_A;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
//Adiabatic cooling
double Adiabatic_Cooling_time(struct temp_ev *pt, struct blob *pt_spec, double R_jet_t){
	return R_jet_t/(pt->v_exp_by_c*vluce_cm);
	//return 1.0/time_adiabatic(pt,pt_spec,pt->R_H_jet_t);
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
//Termine Sistematico + Cooling
double Bfp(double x, struct temp_ev *pt, struct blob *pt_spec){
	return -f_Acc(x,pt) + Cooling(x,pt,pt_spec);
}

//---------------------------------------------------------------------


//---------------------------------------------------------------------
void Wm(double wm, double *WP, double *WM){
	double a,b;
	double abs_wm,ex,ww;
	//printf("Wm wm=%e\n",wm);
	abs_wm=fabs(wm);
	if (abs_wm<0.1){
		a=wm*wm/24;
		b=wm*wm*wm*wm/80;
		ww=1/(1+a+b);
		ex=exp(wm/2);
		*WP=ww*ex;
		*WM=ww/ex;
	}
	else if (wm > 0) {
		ex=exp(wm);
		*WM=wm/(ex-1);
		*WP=*WM+wm;
	}
	else if (wm < 0) {
		ex=exp(-wm);
		*WP=wm/(1-ex);
		*WM=*WP-wm;
	}
	else {
		printf("something wrong in Wm\n");
		exit(0);
	}
	return ;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------



/************************************************************************/
/* funzione che risolve il sistema tri_diagonale                        */
/************************************************************************/
/* il puntatore punta la primo elemento di N_e_1                        */
/************************************************************************/
int solve_sys1(double A[],double B[],double C[],double R[],double u[],unsigned int SIZE){
	//double gam[SIZE];
	double bet,*gam;
	unsigned int j;
	gam=(double*)calloc(SIZE,sizeof(double));

	for(j=1;j<SIZE-1;j++){
		//printf("A[%d]=%e B[%d]=%e C[%d]=%e  \n",j,A[j],j,B[j],j,C[j]);
	}
	

	if(B[0]==0){
		printf("2\n");
		return 2;
	}
	u[0]=R[0]/(bet=B[0]);
	for(j=1;j<SIZE-1;j++){
		gam[j]=C[j-1]/bet;
		bet=B[j]-A[j]*gam[j];
		if(bet==0.0){
			printf("3\n");
			return 3;
		}
		u[j]=(R[j]-A[j]*u[j-1])/bet;
		//printf("1 u[%i]=%e R=%e bet=%e\n",j,u[j],R[j],bet);
	}
	for(j=(SIZE-1);j>=1;j--){
		u[j]-=gam[j+1]*u[j+1];
		//printf("2 u[%i]=%e R=%e \n",j,u[j],R[j]);
	}
	/* test sugli indici */
	for(j=0;j<SIZE;j++){
		//printf("A[%d]=%e B[%d]=%e gam[%d]=%e u=%e \n",j,A[j],j,B[j],j,gam[j],u[j]);
	}
	free(gam);
	//  free_vector(gam,1,pt_griglia_max);
	return 0;
}


