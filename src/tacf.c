/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <R_ext/Applic.h>
#include <math.h>
//#include <stdlib.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
 
#ifndef max
#define max(a,b) ((a < b)?(b):(a))
#endif
#ifndef min
#define min(a,b) ((a < b)?(a):(b))
#endif
 
 

double multip(double *x, double *y, int lx, int ly)  
{
	int i;
	double z = 0;
	if(lx != ly) error("length of x does not match length of y");
	for(i = 0; i < lx; i++) z += x[i]*y[i];
	return z;

}


void tacvfFDWN_C(double *d, int *maxlag, double *x)
{
	int i; 
	double db= (double)(*d);
	x[0] = tgamma(1 - 2*db)/pow(tgamma(1-db), 2);
	for(i = 1; i < (int)(*maxlag)+1; i++) x[i] = ((i - 1 + db)/(i - db))*x[i-1];
}

void tacfFGN_C(double *H, int *maxlag, double *x) {
	int i;
	int ml = (int)(*maxlag);
	double h2 = (double)(*H)*2;
	double temp[ml+2];
	x[0] = 1;
	if(ml > 1) {
		temp[0] = 0;
		temp[1] = 1;
		x[1] = (pow(2, h2) - 2)*0.5;
		for(i = 2; i < ml+2; i++) temp[i] = pow(i, h2);
		for(i = 2; i < ml+1; i++) x[i] = 0.5*(temp[i-1] - 2.0*temp[i] + temp[i+1]);
	}
	else if(ml == 1) x[1] = (pow(2, h2) - 2)*0.5;
}

//From...
double zeta(double s, int n) {
	int k;
	double d[n+1], z, temp; 
	d[0] = 1;
	for(k = 1; k < n+1; k++) {
		temp = tgamma(n+k)/(tgamma(n-k+1)*tgamma(2*k+1));
		d[k] = d[k-1] + n*temp*pow(4.0, k);
	}
	z = 0.0;
	for(k = 0; k < n; k++) {
		z += pow(-1, k)*(d[k] - d[n])/pow(k+1, s);
	}
	z *= (-1)/(d[n]*(1 - pow(2.0, 1 - s)));
	return(z);
}

void tacfHD_C(double *alpha, int *maxlag, double *x) {

	int k, ml = maxlag[0], n = 20; //Hardcode for now.
	double a = alpha[0];
	double c = -1/(2*zeta(a, n));
	x[0] = 1;
	for(k = 1; k < ml + 1; k++) {
		x[k] = c * pow(k, -a);
	}
} 


void tacvfARMA_C(double *phi, int *pp, double *theta, int *qp, int *maxlag, double *res)
{
	int p = (int)*pp;
	int q = (int)*qp;
	int ml = (int)*maxlag;
	int i, j, k, r, info;
	double *C, *b, *theta2, *phi2, *g, *a, *temp;
	int one = 1;

	int *ipiv;
	res[0] = 1;
	for(i = 1; i < (ml+1); i++) res[i] = 0;
	if(max(p, q) > 0) {
		r = max(p, q) +1;
		C = (double *) R_alloc(q+1, sizeof(double));
		b = (double *) R_alloc(r, sizeof(double));
		phi2 = (double *) R_alloc(3*r, sizeof(double));
		theta2 = (double *) R_alloc(q+1, sizeof(double));
		g = (double *) R_alloc(max(r, ml+1), sizeof(double));
		temp = (double *) R_alloc(p, sizeof(double));
		a = (double *) R_alloc(r*r, sizeof(double));
		C[0] = 1; 
		theta2[0] = -1;
		for(i = 0; i < q; i++) theta2[i+1] = theta[i];
		for(i = 0; i < 3*r; i++) phi2[i] = 0;
		phi2[r-1] = -1;
		if(p > 0)
			for(i = 0; i < p; i++) phi2[r+i] = phi[i];
		if(q > 0) {
			for(k = 0; k < q; k++) {
				C[k+1] = -theta[k];
				if(p > 0) for(i = 1; i <= min(p, k+1); i++) C[k+1] += phi[i-1]*C[k+1-i];	
			}
		}
		for(k = 0; k < (q+1); k++) b[k] = 0;
		for(k = 0; k < (q+1); k++)
			for(i = k; i < (q+1); i++)
				b[k] -= theta2[i] * C[i-k];
		for(k = (q+1); k < r; k++) b[k] = 0;

		if(p == 0) {
			for(i = 0; i < (q+1); i++) res[i] = b[i];
			for(i = (q+1); i < (ml+1); i++) res[i] = 0;
		}
		else {
			for(i = 0; i < r; i++) {
				for(j = 0; j < r; j++) {
					if(j == 0) a[i+r*j] = phi2[r+i-1];
					else a[i+r*j] = phi2[r+i-j-1] + phi2[r+i+j-1];
				}
			}

			ipiv = (int *) R_alloc(r, sizeof(int));
			for(i = 0; i < r; i++) b[i] = -b[i];
			F77_CALL(dgesv)(&r, &one, a, &r, ipiv, b, &r, &info);

			if (info < 0)
				error(("argument %d of Lapack routine %s had invalid value"), -info, "dgesv");
			if (info > 0)
				error(("Lapack routine dgesv: system is exactly singular"));

			for(i = 0; i < r; i++) g[i] = b[i];
			if(r <= ml) {
				for(i = r; i < (ml+1); i++) {
					for(j = 1; j <= p; j++) temp[j-1] = g[i-j];
					g[i] = multip(phi, temp, p, p);
				}
			}
			for(i = 0; i < (ml+1); i++) res[i] = g[i];
		}
	}
}


void shift_C(double* x, int *byam, int *n, double* y) 
{
	int i;
	int by = (int)*byam;
	y[0] = x[0];
	for(i = 1; i < *n; i++) y[i*by] = x[i];	
}
