/*****************************************************************************/
/*
/* 16.10.28 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#include "mex.h"
#include "lf4d_matlab.h"
#include "mat_matlab.h"

/* prepare mem before calling lf_4d_resample */
lf_4d_resample_interface(
   double* iLF, /* in: 4D LF */
   int sx,         /* in: origin dim */
   int sy,         /* in: origin dim */
   int nx,
   int ny,
   double* oLF,     /* out: 4D LF */
   int tnx,        /* out: target dim */
   int tny         /* out: target dim */
)
{

float**** iLF4D;
float**** oLF4D;
float** tmp;
int i,j;

/* alloc mem */
ALLOC_4D_LF(&iLF4D,sx,sy,nx,ny);
ALLOC_4D_LF(&oLF4D,sx,sy,tnx,tny);
ALLOC_2D_MATRIX(&tmp,nx,ny);

/*set tmp to zeros*/
/*for(i =0; i<nx; i++)
   for(j =0; j<ny; j++){
      tmp[i][j] = 0.0;
   }
*/
/* copy input */
lf_4d_from_1D(iLF4D,iLF,sx,sy,nx,ny);

/* do resample */
lf_4d_resample(iLF4D, sx,sy,nx,ny,oLF4D,tnx,tny,tmp);
/*copy result back */
lf_4d_to_1D(oLF,oLF4D,sx,sy,tnx,tny);

/* free mem */
FREE_4D_LF(iLF4D,sx,sy,nx,ny);
FREE_4D_LF(oLF4D,sx,sy,tnx,tny);
FREE_2D_MATRIX(tmp,nx,ny);

}


/* The gateway function */
void mexFunction( int noArg, mxArray *poArg[],
                  int niArg, const mxArray *piArg[])
{
    int i;
    int j;

    int nx;
    int ny;
    int sx;
    int sy;

    int tnx;
    int tny;

    double* origin_dim;
    double* target_dim;

    double* iLF;
    double* oLF;


    /* get the value of the scalar input  */
    origin_dim = mxGetPr(piArg[0]);

    sx = (int)origin_dim[0];
    sy = (int)origin_dim[1];
    nx = (int)origin_dim[2];
    ny = (int)origin_dim[3];

    target_dim = mxGetPr(piArg[1]);

    tnx = (int)target_dim[0];
    tny = (int)target_dim[1];

    iLF = mxGetPr(piArg[2]);;

    mexPrintf("Got LF dims: sx %d sy %d nx %d ny %d \n",sx,sy, nx,ny);
    mexPrintf("Will convert to sx%d sy%d nx %d ny %d\n",sx,sy,tnx,tny);

    /* get a pointer to the real data in the output matrix */
    /* create the output matrix */
    mwSize dims[4] = {sy, sx, tny,tnx};
    poArg[0] =  mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);

    oLF = mxGetPr(poArg[0]);

    lf_4d_resample_interface(iLF,sx,sy,nx,ny,oLF,tnx,tny);
}
