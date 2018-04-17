/*****************************************************************************/
/*
/* 16.10.28 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#include "mex.h"
#include "lf4d_matlab.h"
#include "mat_matlab.h"



img_2D_from_1D(float** iIMG2D,double* iIMG,int nx,int ny,int bx,int by){
   int i,j;
   int H = ny+2*by;

   for(i = bx; i<nx+bx; i++)
      for(j =by; j<ny+by; j++){
         iIMG2D[i][j] = (float)iIMG[(i-bx)*ny + j-by];
      }

}

img_2D_to_1D(double* oIMG,float** oIMG2D,int nx,int ny,int bx,int by){
   int i,j;
   int H = ny+2*by;

   for(i = bx; i<nx+bx; i++)
      for(j =by; j<ny+by; j++){
         oIMG[(i-bx)*ny +(j-by)] = oIMG2D[i][j];
      }

}


/* prepare mem before calling lf_4d_resample */
img_resample_interface(
   double* iIMG, /* in: 4D LF */
   int nx,
   int ny,
   double* oIMG,     /* out: 4D LF */
   int tnx,        /* out: target dim */
   int tny         /* out: target dim */
)
{

float** iIMG2D;
float** oIMG2D;
float** tmp;
int i,j;

int tmpnx,tmpny;

tmpnx = (tnx<nx)?nx:tnx;
tmpny = (tny<ny)?ny:tny;

/* alloc mem */
ALLOC_2D_MATRIX(&iIMG2D,nx,ny);
ALLOC_2D_MATRIX(&oIMG2D,tnx,tny);
ALLOC_2D_MATRIX(&tmp,tmpnx,tmpny);

/*copy input*/
img_2D_from_1D(iIMG2D,iIMG,nx,ny,0,0);

/* set tmp to zeros */
for(i=0; i<nx; i++)
   for(j =0; j<ny; j++){
      tmp[i][j] = 0.0;
   }

for(i=0; i<tnx; i++)
   for(j =0; j<tny; j++){
      oIMG2D[i][j] = 0.0;
   }

/* do resample */

/*my_resample_2d_y(iIMG2D, nx,ny,bx,by,oIMG2D,tny);*/
my_resample_2d(iIMG2D, nx,ny,0,0,oIMG2D,tnx,tny,tmp);
/*lf_4d_resample(iLF4D, nx,sy,oIMG2D,tnx,tny,tmp);*/

/*copy result back */
img_2D_to_1D(oIMG,oIMG2D,tnx,tny,0,0);

/* free mem */
FREE_2D_MATRIX(tmp,tmpnx,tmpny);
FREE_2D_MATRIX(oIMG2D,tnx,tny);
FREE_2D_MATRIX(iIMG2D,nx,ny);

}


/* The gateway function */
void mexFunction( int noArg, mxArray *poArg[],
                  int niArg, const mxArray *piArg[])
{
    int i;
    int j;

    int nx;
    int ny;

    int tnx;
    int tny;

    double* origin_dim;
    double* target_dim;

    double* iIMG;
    double* oIMG;


    /* get the value of the scalar input  */
    origin_dim = mxGetPr(piArg[0]);

    nx = (int)origin_dim[0];
    ny = (int)origin_dim[1];

    target_dim = mxGetPr(piArg[1]);

    tnx = (int)target_dim[0];
    tny = (int)target_dim[1];

    iIMG = mxGetPr(piArg[2]);;

    mexPrintf("Got IMG dims: nx %d ny %d \n",nx,ny);
    mexPrintf("Will convert to nx %d ny %d\n",tnx,tny);

    /* get a pointer to the real data in the output matrix */
    /* create the output matrix */
    mwSize dims[2] = {tny,tnx};
    poArg[0] =  mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    oIMG = mxGetPr(poArg[0]);

    img_resample_interface(iIMG,nx,ny,oIMG,tnx,tny);
}
