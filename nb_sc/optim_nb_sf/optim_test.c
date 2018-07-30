#include <stdio.h>
#include <assert.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

void lbfgsb_(int n, int m, double *x, double *l, double *u, int *nbd,
       double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
       void *ex, double factr, double pgtol,
       int *fncount, int *grcount, int maxit, char *msg,
       int trace, int nREPORT);

const int n = 2;
const int m = 2;

struct nb_dataT {
   int n;
   double *k;
   double *sf;
};

struct nb_dataT test_nb_data =
    ( struct nb_dataT ){
      10,
      (double[]){ 0, 0, 0, 0, 0, 0, 1, 2, 6, 6 },
      (double[]){ 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.30, 1.50, 1.90, 2.00 } };

double nll_nb_sf( int n, double *x, void *p ) {
   /* mu in x[0], disp in x[1] */
   struct nb_dataT * data = (struct nb_dataT *) p;
   double sum = 0;
   for(int i = 0; i < data->n; ++i) {
     sum += lgammafn( data->k[i] + 1. / x[1] ) - lgammafn( 1. / x[1] ) - lgammafn( data->k[i] + 1. ) - 
        ( data->k[i] + 1. / x[1]) * log( 1 + data->sf[i] * x[1] * x[0] ) + data->k[i] * log( data->sf[i] * x[1] * x[0] );
   }
   printf( "%f,%f -> %f\n", x[0], x[1], sum );
   return( -sum );
}

void gr_nll_nb_sf( int n, double *x, double *g, void *p ) {
   struct nb_dataT * data = (struct nb_dataT *) p;
   double sum0 = 0;
   double sum1 = 0;
   int i;
   for( i = 0; i < data->n; i++ ) {
      sum0 += ( data->k[i] - data->sf[i] * x[0] ) / ( x[0] + data->sf[i] * x[1] * x[0]*x[0] );
      sum1 += ( x[1] * ( data->k[i] - data->sf[i] * x[0] ) / ( 1. + data->sf[i] * x[1] * x[0] ) + 
         log( 1. + data->sf[i] * x[1] * x[0] ) - digamma( data->k[i] + 1/x[1] ) + digamma( 1./x[1] ) ) / (x[1]*x[1]); }      
   g[0] = -sum0;
   g[1] = -sum1;
}

int main() {
   double x[] = { .5, 1 };
   int nbd[] = { 1, 1 };
   double l[] = { 0.1, 0.1 };
   double fmin, g[2];
   int fail, fncount, grcount;
   char msg[200];

   lbfgsb_( 
      n, /* int n */ 
      m, /* int m */ 
      x, /* double *x */
      l, /* double *l */
      NULL, /* double *u */
      nbd, /* int *nbd */
      &fmin, /* double *Fmin */ 
      &nll_nb_sf, /* optimfn fn */
      &gr_nll_nb_sf, /* optimgr gr */
      &fail, /* int *fail */
      &test_nb_data, /* void *ex */
      1e10, /* double factr */
      1e-10, /* double pgtol */
      &fncount, /* int *fncount */
      &grcount, /* int *grcount */
      100, /* int maxit */ 
      msg, /* char *msg */
      0, /* int trace */
      1 /* int nREPORT */ );

      if( fail )
         printf( "failed!\n");

      printf( "Optimum at %f, %f\n", x[0], x[1] );
      printf( "function value %f\n", fmin );
}