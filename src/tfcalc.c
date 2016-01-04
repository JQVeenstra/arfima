#include <math.h>

#define max(a,b) a > b ? a : b
#define min(a,b) a > b ? b : a



void tfcalc(double *y, int *n, double *x, double *delta, int *r, double *omega, int *s, int *b, int *num, double *meanval) 
{

	int i, j, u, k;
	int nn = num[0], n0 = n[0];
	double yh[n0*nn]; 
	int rr = 0, ss = 0;
	double mean = meanval[0];
	for(i = 0; i < n0*nn; i++) yh[i] = -mean; //*******************
	for(k = 0; k < nn; k++) {

		u = max(r[k], s[k]+b[k]) - 1;
		
		for(i = u; i < n0; i++) {
			
			for(j = 0; j < r[k]; j++) {
				
				yh[i+k*n0] += delta[rr+j]*yh[i - j - 1 + k*n0]; 
			
			}
			
			yh[i + k*n0] += omega[ss+0]*x[i-b[k] + k*n0];
			
			for(j = 1; j < s[k]; j++) {
			
				yh[i + k*n0] -= omega[ss+j]*x[i - j - b[k] + k*n0];
			
			}
		
		}
		
		rr += r[k];
		ss += s[k];
	}
	for(k = 0; k < nn; k++) {
		for(i = 0; i < n0; i++) {
			
			y[i] -= yh[i+k*n0];
			
		}	
	}
}

