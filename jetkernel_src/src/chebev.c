#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

double chebev(double a, double b, double c[], int m, double x)
{

	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0){
		exit(0);
		printf("NRERROR x not in range in routine chebev");
	}
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}
/* (C) Copr. 1986-92 Numerical Recipes Software !'K4$<%#110L&")|oV'4. */
