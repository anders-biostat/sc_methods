#include <stdio.h>
#include <stdlib.h>

#include <R_ext/Applic.h>
#include <R_ext/Print.h>	/* for Rprintf */
#include <R_ext/Arith.h>
#include <R_ext/Memory.h>

#include "lbfgsb.c"

char * _( char * s ) {
   return s;
}

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}

void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
	    void *ex, double factr, double pgtol,
	    int *fncount, int *grcount, int maxit, char *msg,
	    int trace, int nREPORT)
{
    char task[60];
    double f, *g,  *wa;
    int tr = -1, iter = 0, *iwa, isave[21];
    isave[12] = 0; // -Wall

    if(n == 0) { /* not handled in setulb */
	*fncount = 1;
	*grcount = 0;
	*Fmin = fminfn(n, u, ex);
	strcpy(msg, "NOTHING TO DO");
	*fail = 0;
	return;
    }
    if (nREPORT <= 0)
	error(_("REPORT must be > 0 (method = \"L-BFGS-B\")"));
    switch(trace) {
    case 2: tr = 0; break;
    case 3: tr = nREPORT; break;
    case 4: tr = 99; break;
    case 5: tr = 100; break;
    case 6: tr = 101; break;
    default: tr = -1; break;
    }

    *fail = 0;
    g = vect(n);
    /* this needs to be zeroed for snd in mainlb to be zeroed */
    wa = (double *) S_alloc(2*m*n+4*n+11*m*m+8*m, sizeof(double));
    iwa = (int *) R_alloc(3*n, sizeof(int));
    strcpy(task, "START");
    while(1) {
	setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
	       tr, isave);
/*	Rprintf("in lbfgsb - %s\n", task);*/
	if (strncmp(task, "FG", 2) == 0) {
	    f = fminfn(n, x, ex);
	    if (!R_FINITE(f))
		error(_("L-BFGS-B needs finite values of 'fn'"));
	    fmingr(n, x, g, ex);
	} else if (strncmp(task, "NEW_X", 5) == 0) {
	    iter++;
	    if(trace == 1 && (iter % nREPORT == 0)) {
		Rprintf("iter %4d value %f\n", iter, f);
	    }
	    if (iter > maxit) {
		*fail = 1;
		break;
	    }
	} else if (strncmp(task, "WARN", 4) == 0) {
	    *fail = 51;
	    break;
	} else if (strncmp(task, "CONV", 4) == 0) {
	    break;
	} else if (strncmp(task, "ERROR", 5) == 0) {
	    *fail = 52;
	    break;
	} else { /* some other condition that is not supposed to happen */
	    *fail = 52;
	    break;
	}
    }
    *Fmin = f;
    *fncount = *grcount = isave[12];
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit && *fail == 0) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    strcpy(msg, task);
}
