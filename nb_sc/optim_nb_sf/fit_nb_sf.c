#include <stdio.h>
#include <assert.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

struct nb_dataT {
   int n;
   double *k;
   double *sf;
};

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

SEXP fit_nb_sf( SEXP k, SEXP sf ) {
   double x[] = { .5, 1 };  // initial values
   int nbd[] = { 1, 1 };    // lower bounds only
   double l[] = { 1e-15, 1e-15 };  // lower bounds
   double fmin, g[2];
   int fail, fncount, grcount;
   char msg[200];

   if( Rf_length(k) != Rf_length(sf) )
   	  Rf_error( " Lengths are different.");

   struct nb_dataT data;
   data.n = Rf_length(k);
   data.k = REAL(k);
   data.sf = REAL(sf);

   lbfgsb( 
      2, /* int n */ 
      2, /* int m */ 
      x, /* double *x */
      l, /* double *l */
      NULL, /* double *u */
      nbd, /* int *nbd */
      &fmin, /* double *Fmin */ 
      &nll_nb_sf, /* optimfn fn */
      &gr_nll_nb_sf, /* optimgr gr */
      &fail, /* int *fail */
      &data, /* void *ex */
      1e10, /* double factr */
      1e-10, /* double pgtol */
      &fncount, /* int *fncount */
      &grcount, /* int *grcount */
      100, /* int maxit */ 
      msg, /* char *msg */
      0, /* int trace */
      1 /* int nREPORT */ );

   if( fail )
      Rprintf( "failed!\n");

   Rprintf( "Optimum at %f, %f\n", x[0], x[1] );
   Rprintf( "function value %f\n", fmin );

   SEXP res = Rf_allocVector( REALSXP, 2 );
   REAL(res)[0] = x[0];
   REAL(res)[1] = x[1];
   return res;
}

R_CallMethodDef callMethods[] = {
   { "fit_nb_sf", (DL_FUNC) &fit_nb_sf, 2 },
   { NULL, NULL, 0 }
};

void R_init_optim_test( DllInfo *info )
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
