#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Applic.h>

void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *f, double *g, double factr, double *pgtol,
	    double *wa, int * iwa, char *task, int iprint, int *isave);

void lbfgsb_(int n, int m, double *x, double *l, double *u, int *nbd,
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
	printf("REPORT must be > 0 (method = \"L-BFGS-B\")");
    switch(trace) {
    case 2: tr = 0; break;
    case 3: tr = nREPORT; break;
    case 4: tr = 99; break;
    case 5: tr = 100; break;
    case 6: tr = 101; break;
    default: tr = -1; break;
    }

    *fail = 0;
    g = (double *) calloc ( n, sizeof(double) );
    /* this needs to be zeroed for snd in mainlb to be zeroed */
    wa = (double *) calloc(2*m*n+4*n+11*m*m+8*m, sizeof(double));
    iwa = (int *) calloc(3*n, sizeof(int));
    strcpy(task, "START");
    while(1) {
	setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
	       tr, isave);
	/*printf("in lbfgsb - %s\n", task);*/
	if (strncmp(task, "FG", 2) == 0) {
	    f = fminfn(n, x, ex);
	    fmingr(n, x, g, ex);
	} else if (strncmp(task, "NEW_X", 5) == 0) {
	    iter++;
	    if(trace == 1 && (iter % nREPORT == 0)) {
		printf("iter %4d value %f\n", iter, f);
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
	printf("final  value %f \n", *Fmin);
	if (iter < maxit && *fail == 0) printf("converged\n");
	else printf("stopped after %i iterations\n", iter);
    }
    strcpy(msg, task);
}

