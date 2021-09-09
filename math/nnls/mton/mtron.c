

/*
 * see mtron.m for helptext
 */


#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef NULL
#define NULL   0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef MAX
#define MAX(a,b)         ((a < b) ?  (b) : (a))
#endif

#ifndef MIN
#define MIN(a,b)         ((a > b) ?  (b) : (a))
#endif

#ifndef ABS
#define ABS(x)           (((x) < 0) ? -(x) : (x))
#endif



/*
 * additional storage for preconditioner.
 */
#define TRON_NNZPLUS_FAC   30


/*
 * Prototype for TRON routine.
 */
void dtron_(int *n, double *x, double *xl, double *xu, 
        double *f, double *G, double *H, double *Hdiag, int *jc, int *ir, 
	    double *frtol, double *fatol, double *fmin, double *cgtol, 
        int *itermax, double *delta, char *task,
	    double *b, double *bdiag, int *bcol_ptr, int *brow_ind,
	    double *l, double *ldiag, int *lcol_ptr, int *lrow_ind,
	    double *xc, double *s, int *indfree, int *isave, double *dsave, 
	    double *wa, int *iwa); 


/*
 * Helper routine to copy f
 */ 
void copy_f(double *f, const mxArray *mxf)
{
  double *ftemp;
  
  /* check whether mxf corresponds to a scalar double,
       and whether the data can be extracted. */
  ftemp = mxGetPr(mxf);
  if ((ftemp == NULL) || (mxGetM(mxf) * mxGetN(mxf) != 1))
    mexErrMsgTxt("mtron: error in copy_f, mxf must be scalar.");
  
  /* copy data. */
  f[0] = ftemp[0];
}

/*
 * Helper Routine to copy G
 */
void copy_G(int n, double *G, const mxArray *mxG)
{
    double *Gtemp;
    int i;
    
    /* check whether the dimensions match */
    if ( (mxGetM(mxG) != n) || (mxGetN(mxG) != 1) )
        mexErrMsgTxt("mtron: dimension mismatch in copy_G");
    
    /* extract data and check for correctness */
    Gtemp = mxGetPr(mxG);
    if (Gtemp == NULL)
        mexErrMsgTxt("mtron: NULL pointer Gtemp in copy_G");
    
    /* copy data. */
    for (i = 0; i < n; i++)
        G[i] = Gtemp[i];
}


/*
 * Helper routine to copy H and Hdiag
 */
void copy_H(int n, int nnz, 
            double *H, double *Hdiag, int *ir, int *jc, 
            const mxArray *mxH, const mxArray *mxHdiag)
{
    double *Ht, *Hdiagt;
    int *irt, *jct;
    int i;
    
    /* first check dimensions. */
    if ( (mxGetM(mxH) != n) || (mxGetN(mxH) != n) )
        mexErrMsgTxt("mtron: dimension mismatch for mxH in copy_H");
    
    /* extract ir, jc and H temporary arrays. */
    irt = mxGetIr(mxH);
    jct = mxGetJc(mxH);
    Ht = mxGetPr(mxH);
    /* check whether pointers are valid */
    if ( (irt == NULL) || (jct == NULL) || (Ht == NULL) )
        mexErrMsgTxt("mtron: NULL pointer in copy_H (1)");
    
    /* check whether nnz is an upper bound for the number of elements! */
    if (jct[n] > nnz)
        mexErrMsgTxt("mtron: in copy_H jc[end] > nnz!!!");

    /* copy ir, jc, H
       indices are shifted by one (fortran to C conversion) */
    for (i = 0; i < n+1; i++)
        jc[i] = jct[i] + 1;
    for (i = 0; i < jc[n]; i++) {
        ir[i] = irt[i] + 1;
        H[i] = Ht[i];
    }

    /* check dimensions for Hdiag */
    if ( (mxGetM(mxHdiag) != n) || (mxGetN(mxHdiag) != 1) )
        mexErrMsgTxt("mtron: dimension mismatch for mxHdiag in copy_H");
    /* copy Hdiag */
    Hdiagt = mxGetPr(mxHdiag);
    if (Hdiagt == NULL)
        mexErrMsgTxt("mtron: NULL pointer in copy_H (2)");
    /* copy Hdiag */
    for (i = 0; i < n; i++)
        Hdiag[i] = Hdiagt[i];
    
}


void create_info(mxArray **info, double delta, int ncgit)
{
  mxArray *temp;
  
  static const int dims[] = {1};
  static const char *names[] = { "delta", "ncgit" };
  
  /* create the structure */
  *info = mxCreateStructArray(1, dims, 2, names);
  if (*info == NULL) 
    mexErrMsgTxt("mtron:create_info: Couldn't create info field.");
  
  /* create the delta field */
  temp = mxCreateDoubleMatrix(1, 1, mxREAL);
  *(mxGetPr(temp)) = delta;
  mxSetFieldByNumber(*info, 0, 0, temp);
  
  /* create the ncgit field */
  temp = mxCreateDoubleMatrix(1, 1, mxREAL);
  *(mxGetPr(temp)) = ((double)(ncgit));
  mxSetFieldByNumber(*info, 0, 1, temp);
  
}


/*
 * Actual MEX routine.
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{

    /*****************************************************************/
    /*                      LOCAL VARIABLES                          */
    /*****************************************************************/
    
    char msg[20];              /* input message from Matlab call        */
    int i;                     /* just a counting index                 */
    double *xtemp;             /* temporary variable for copying x.     */
     
    static char task[100];     /* TRON message passing                  */
    
    static int isInitialized,  /* initialization flag                   */
               nnz,            /* max. number of non-zeros in the hessian*/
               max_nnz,        /* max. nnz in L matrix                  */
               itermax,        /* max. num. CG iterations               */
               n;              /* problem dimension                     */
                
    static double cgtol,          /* termination tolerance for CG       */
                  fmin,           /* lower bound for objective function */
                  delta,          /* trust region radius                */
                  frtol, fatol,   /* relative and absolute tolerance    */
                  *xl, *xu;       /* lower and upper bound for x        */
                  
    
    /* definition of the TRON working arrays (see TRON help text)
     * f : current function value
     * G : current gradient
     * H, Hdiag : current Hessian
     * b ... (n x n), length nnz
     * l ... (n x n), length nnz + n * p = max_nnz (see below)
     * xc : dim n, s: n, indfree: n, isave : 3, wa : 7 * n, iwa : 3 * n */
    static double *x, *f, *G, *H, *Hdiag,
                  *b, *bdiag, *l, *ldiag,
                  *xc, *s, *dsave, *wa;
    
    static int *ir, *jc, *bcol_ptr, *brow_ind, *lcol_ptr, *lrow_ind, 
               *indfree, *isave, *iwa;
    
    
    /*****************************************************************/
    /*                  READ COMMAND STRING                          */
    /*****************************************************************/
    
    
    /* First input check: what type of input? */
    
    /* if the number of input arguments is 0 then return an error. */
    if (nrhs < 1) 
        mexErrMsgTxt("mtron: Must have at least one input argument!");
    /* if the first argument is not a string, return an error. */
    if (!mxIsChar(prhs[0]))
        mexErrMsgTxt("mtron: First argument must be a string.");
    /* first argument must have length at least 1  */
    if (mxGetM(prhs[0]) * mxGetN(prhs[0]) == 0) 
        mexErrMsgTxt("mtron: First argument must be a string of length at least one.");
    
    /* if at least one input argument is present, read it out .
     * If the message is longer than 20 characters, it doesnt
     * matter, since we use only the first one.  */
    mxGetString(prhs[0], msg, 20);
    
    
    /*****************************************************************/
    /*                    CLEANUP COMMAND                            */
    /*****************************************************************/
    /* this is called if msg = 'CLEAN' or if msg = 'INIT' and
     * the the temporary arrays are still allocated... */
    if ( (isInitialized) && ((msg[0] == 'C') || (msg[0] == 'I')) ) {
        
        /* clear all memory persistent variables. */
        mxFree(x);
        mxFree(xl);
        mxFree(xu);
        
        mxFree(f);
        mxFree(G);
        mxFree(H);
        mxFree(Hdiag);
        mxFree(ir);
        mxFree(jc);
        
        mxFree(b);
        mxFree(bdiag);
        mxFree(bcol_ptr);
        mxFree(brow_ind);
        mxFree(l);
        mxFree(ldiag);
        mxFree(lcol_ptr);
        mxFree(lrow_ind);
        mxFree(xc);
        mxFree(s);
        mxFree(indfree);
        mxFree(isave);
        mxFree(dsave);
        mxFree(wa);
        mxFree(iwa);
        
        /* reset the initialization flag */
        isInitialized = FALSE;
    }
    
    
    /*****************************************************************/
    /*                      INIT COMMAND                             */
    /*****************************************************************/
    if (msg[0] == 'I') {

        /* check input: need f, G, H, Hdiag. */
        if (nrhs < 13) 
            mexErrMsgTxt("mtron(INIT): Invalid number of arguments. See 'help mtron'.");
        
        /* read all parameters  */
        /* problem dimension = length of x0 array  */
        n = mxGetM(prhs[1]);
        /* fmin = lower bound for F  */
        fmin = *( mxGetPr(prhs[4]) );
        /* delta = initial trust region radius  */
        delta = *( mxGetPr(prhs[9]) );
        /* frtol, fatol, cgtol : termination tolerance values  */
        frtol = *( mxGetPr(prhs[10]) );
        fatol = *( mxGetPr(prhs[11]) );
        cgtol = *( mxGetPr(prhs[12]) );
        /* nnz = maximum number of non-zeros in the hessian matrix
         * prhs[7] is the hessian. get the jc array and extract nnz
         * from it. Just in case, we should maybe take a bit more ??TODO?*/
        jc = mxGetJc(prhs[7]);
        nnz = (jc[n] + 100) * 2;
        /* max_nnz is computed by the usual formula. */
        max_nnz = nnz + TRON_NNZPLUS_FAC * n;
        /* set max. CG iterations to default value ***TODO*** */
        itermax = n;
      
        /* allocate the working arrays and make them memory persistent */
        x = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(x);
        xl = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(xl);
        xu = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(xu);
        
        f = mxCalloc(1, sizeof(double)); 
        mexMakeMemoryPersistent(f);
        G = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(G);
        H = mxCalloc(nnz, sizeof(double));
        mexMakeMemoryPersistent(H);
        Hdiag = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(Hdiag);
        ir = mxCalloc(nnz, sizeof(int));
        mexMakeMemoryPersistent(ir);
        jc = mxCalloc(n+1, sizeof(int));
        mexMakeMemoryPersistent(jc);
        
        b = mxCalloc(nnz, sizeof(double));
        mexMakeMemoryPersistent(b);
        bdiag = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(bdiag);
        bcol_ptr = mxCalloc((n+1), sizeof(int));
        mexMakeMemoryPersistent(bcol_ptr);
        brow_ind = mxCalloc(nnz, sizeof(int));
        mexMakeMemoryPersistent(brow_ind);
        l = mxCalloc(max_nnz, sizeof(double));
        mexMakeMemoryPersistent(l);
        ldiag = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(ldiag);
        lcol_ptr = mxCalloc((n+1), sizeof(int));
        mexMakeMemoryPersistent(lcol_ptr);
        lrow_ind = mxCalloc(max_nnz, sizeof(int));
        mexMakeMemoryPersistent(lrow_ind);
        xc = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(xc);
        s = mxCalloc(n, sizeof(double));
        mexMakeMemoryPersistent(s);
        indfree = mxCalloc(n, sizeof(int));
        mexMakeMemoryPersistent(indfree);
        isave = mxCalloc(3, sizeof(int));
        mexMakeMemoryPersistent(isave);
        dsave = mxCalloc(3, sizeof(double));
        mexMakeMemoryPersistent(dsave);
        wa = mxCalloc(7 * n, sizeof(double));
        mexMakeMemoryPersistent(wa);
        iwa = mxCalloc(3 * n, sizeof(int));
        mexMakeMemoryPersistent(iwa);
        
        /* set initialized flag */
        isInitialized = TRUE;
        
        /* copy x, xl, xu
         * (here, I abuse the function copy_G)  */
        copy_G(n, x, prhs[1]);
        copy_G(n, xl, prhs[2]);
        copy_G(n, xu, prhs[3]);
        
        
        /* copy problem information into f, G, H, Hdiag  */
        copy_f(f, prhs[5]);
        copy_G(n, G, prhs[6]);
        copy_H(n, nnz, H, Hdiag, ir, jc, prhs[7], prhs[8]);
        
        /* set task to 'START'  */
        sprintf(task, "START");
        
        /* call TRON.  */
        dtron_(&n, x, xl, xu, f, G, H, Hdiag, jc, ir,
               &frtol, &fatol, &fmin, &cgtol, &itermax, &delta, task,
               b, bdiag, bcol_ptr, brow_ind,
               l, ldiag, lcol_ptr, lrow_ind,
               xc, s, indfree, isave, dsave, wa, iwa);

        /* write output  */

        /* copy the task string */
        if (nlhs > 0)
            plhs[0] = mxCreateString(task);
        
        /* copy the current iterate */
        if (nlhs > 1) {
            plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
            xtemp = mxGetPr(plhs[1]);
            for (i = 0; i < n; i++)
                xtemp[i] = x[i];
        }

        /* create an info field  */
        if (nlhs > 2)
          create_info(&(plhs[2]), delta, isave[2]);
/*         {   plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        *(mxGetPr(plhs[2])) = delta; } */
        
    }
    
    
    /*****************************************************************/
    /*                      STEP COMMAND                             */
    /*****************************************************************/
    if (msg[0] == 'S') {
        
        if (isInitialized != 1)
            mexErrMsgTxt("mtron(STEP): interface has not been initialized. Call mtron('INIT') first!");
        
        /* copy problem information into f, G, H, Hdiag */
        switch (task[0]) {
            case 'F':
                if (nrhs != 2)
                    mexErrMsgTxt("mtron(STEP): task = 'F' >>> need 2 parameters!");
                copy_f(f, prhs[1]);
                break;
                
            case 'G':
                if (nrhs != 4)
                    mexErrMsgTxt("mtron(STEP): task = 'GH' >>> need 4 parameters!");
                copy_G(n, G, prhs[1]);
                copy_H(n, nnz, H, Hdiag, ir, jc, prhs[2], prhs[3]);
                break;
                
            case 'N':
                /* do nothing if x was updated, but still check input */
                if (nrhs != 1)
                    mexErrMsgTxt("mtron(STEP): task = 'NEWX' >>> need 1 parameter!");
                break;
                
            case 'C':
                /* mtron was called even though the last message was
                 * that it has already converged.  */
                mexErrMsgTxt("mtron(STEP): TRON has already converged!");
                
            case 'W':
                /* mtron was called even though the last message was a
                 * a warning, i.e., TRON has terminated without succes */
                mexErrMsgTxt("mtron(STEP): TRON has already terminated with a warning!");
            
        }
        
        /* call TRON */
        dtron_(&n, x, xl, xu, f, G, H, Hdiag, jc, ir,
               &frtol, &fatol, &fmin, &cgtol, &itermax, &delta, task,
               b, bdiag, bcol_ptr, brow_ind,
               l, ldiag, lcol_ptr, lrow_ind,
               xc, s, indfree, isave, dsave, wa, iwa);
        
        /* copy the task string */
        if (nlhs > 0)
            plhs[0] = mxCreateString(task);
        
        /* copy the current iterate */
        if (nlhs > 1) {
            plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
            xtemp = mxGetPr(plhs[1]);
            for (i = 0; i < n; i++)
                xtemp[i] = x[i];
        }
        
        /* create an info */
        if (nlhs > 2) 
          create_info(&(plhs[2]), delta, isave[2]);
/*        {
            plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        *(mxGetPr(plhs[2])) = delta; } */
        
    }
    
    /*****************************************************************/
    /*                 FREE (INDFREE) OPTION                         */
    /*****************************************************************/
    if (msg[0] == 'F') {
      if (isInitialized != 1)
        mexErrMsgTxt("mtron(FREE): interface has not been initialized. Call mtron('INIT') first!");
      /* check number of input/output variables */
      if ( (nrhs != 1) || (nlhs != 1) )
        mexErrMsgTxt("mtron(FREE): requires 1 input, 1 output argument.");
      /* create output array. */
      plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
      xtemp = mxGetPr(plhs[0]);
      /* copy data */
      for (i = 0; i < n; i++)
        xtemp[i] = (double)(indfree[i]);
      
    }

    
    /*****************************************************************/
    /*                 UNRECOGNIZED OPTION                           */
    /*****************************************************************/
    if ((msg[0] != 'I') && (msg[0] != 'C') && (msg[0] != 'S')
        && (msg[0] != 'F') ) {
        mexErrMsgTxt("mtron: Unrecognized command option");
    }
    
}

