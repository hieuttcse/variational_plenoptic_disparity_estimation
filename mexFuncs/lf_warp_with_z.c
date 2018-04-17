/*****************************************************************************/
/*
/* 16.10.28 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#include "mex.h"
#include "lf4d_matlab.h"
#include "mat_matlab.h"




void set_bounds_2d

(
                     /********************************************************/
    float **A,       /* image matrix                                         */
    int   nx,        /* size in x direction                                  */
    int   ny,        /* size in y direction                                  */
    int   bx,        /* boundary in x direction                              */
    int   by,        /* boundary in y direction                              */
    float a          /* set boundaries to a                                  */
                     /********************************************************/
)

/* set boundaries of the 2-D array A to value a */

{
                     /********************************************************/
int  i,j ;           /* loop variables                                       */
                     /********************************************************/

/* upper and lower boundary */
 for (i=0; i<nx+2*bx; i++)
     for (j=1; j<=by; j++)
     {
	 A[i][by   -j]=a;
	 A[i][ny+by-1+j]=a;
     }

/* left and right boundary */
 for (i=1; i<=bx; i++)
     for (j=by; j<ny+by; j++)
     {
	 A[bx     -i][j]=a;
	 A[nx+bx-1+i][j]=a;
     }

return;
}

/*---------------------------------------------------------------------------*/


void mirror_bounds_2d
(
                     /********************************************************/
    float **A,       /* image matrix                                         */
    int   nx,        /* size in x direction                                  */
    int   ny,        /* size in y direction                                  */
    int   bx,        /* boundary in x direction                              */
    int   by         /* boundary in y direction                              */
                     /********************************************************/
)

/* mirror boundaries of the 2-D array A */

{
                     /********************************************************/
int  i,j ;           /* loop variables                                       */
                     /********************************************************/


/* upper and lower boundary */
 for (i=bx; i<nx+bx; i++)
     for (j=1; j<=by; j++)
     {
	 A[i][by     -j]=A[i][by-1 +j];
	 A[i][ny+by-1+j]=A[i][ny+by-j];
     }

/* left and right boundary */
 for (i=1; i<=bx; i++)
     for (j=0; j<ny+2*by; j++)
     {
	 A[bx     -i][j]=A[bx-1 +i][j];
	 A[nx+bx-1+i][j]=A[nx+bx-i][j];
     }

return;
}


void view_backward_registration
(
        /**************************************************/
        float **f1,        /* in  : center image                                */
        float **f2,        /* in  : sub-apeture image                                */
        float **f2_bw,     /* out : sub-apeture image (motion compensated)           */
        float **z,         /* in     : x-component of displacement field     */
        int s,             /* in: view index s*/
        int t,             /* in: view index t*/
        int nx,            /* in     : size in x-direction                   */
        int ny,            /* in     : size in y-direction                   */
        int bx,            /* in     : boundary size in x-direction          */
        int by,            /* in     : boundary size in y-direction          */
        float hx,          /* in     : grid spacing in x-direction           */
        float hy           /* in     : grid spacing in y-direction           */
                           /**************************************************/
)

/* creates warped version of image f2 by means of bilinear interpolation */

{
        /**************************************************/
        int i,j;           /* loop variables                                 */
        int ii,jj;         /* pixel coordinates                              */
        float ii_fp,jj_fp; /* subpixel coordinates                           */
        float delta_i,delta_j; /* subpixel displacement                          */
        float hx_1,hy_1;   /* time saver                                     */
                           /**************************************************/

/* compute time savers */
        hx_1=1.0/hx;
        hy_1=1.0/hy;

/* set boundaries zero */
        set_bounds_2d(f2,nx,ny,bx,by,(float)0.0);

        for (i=bx; i<nx+bx; i++)
                for (j=by; j<ny+by; j++)
                {
                        /* compute subpixel location */
                        ii_fp=i+(z[i][j]*s*hx_1);
                        jj_fp=j+(z[i][j]*t*hy_1);

                        /* if the required image information is out of bounds */
                        if ((ii_fp<bx)||(jj_fp<by)||(ii_fp>(nx+bx-1))||(jj_fp>(ny+by-1)))
                        {
                                /* assume zero flow, i.e. set warped 2nd image to 1st image */
                                f2_bw[i][j]=f1[i][j];

                        }
                        /* if required image information is available */
                        else
                        {
                                /* compute integer index of upper left pixel */
                                ii=(int)floor(ii_fp);
                                jj=(int)floor(jj_fp);

                                /* compute subpixel displacement */
                                delta_i = ii_fp-(float)ii;
                                delta_j = jj_fp-(float)jj;

                                /* perform bilinear interpolation */
                                f2_bw[i][j]   = (1.0-delta_i)*(1.0-delta_j) * f2[ii  ][jj  ]
                                                +      delta_i *(1.0-delta_j) * f2[ii+1][jj  ]
                                                + (1.0-delta_i)*     delta_j  * f2[ii  ][jj+1]
                                                +      delta_i *     delta_j  * f2[ii+1][jj+1];
                        }
                }
}



void lf_backward_registration(
   float **** lf_res,
   float **** lf_res_warp,
   float ** z,
   int sx,
   int sy,
   int nx_fine,
   int ny_fine,
   int bx,
   int by,
   float hx_fine,
   float hy_fine
)
{
   int cs,ct; /*view center s,t*/
   int i,j,s,t;
   float **f0;
   float **f1;
   float ** f1_warp;


  ALLOC_2D_MATRIX(&f0,nx_fine+2*bx,ny_fine+2*by);
  ALLOC_2D_MATRIX(&f1,nx_fine+2*bx,ny_fine+2*by);
  ALLOC_2D_MATRIX(&f1_warp,nx_fine+2*bx,ny_fine+2*by);


   ct = (int) floor((sy-1)/2);
   cs = (int) floor((sx-1)/2);


  lf_4d_extractview_with_boundary(f0,lf_res,sx,sy,nx_fine,ny_fine,bx,by,cs,ct);
  /*mirror_bounds_2d(f0,nx,ny,bx,by);*/
  mirror_bounds_2d(f0,nx_fine,ny_fine,bx,by);
  mirror_bounds_2d(z,nx_fine,ny_fine,bx,by);


   for( s=-cs; s<=cs; s++)
      for(t=-ct; t<=ct; t++){
         lf_4d_extractview_with_boundary(f1,lf_res,sx,sy,nx_fine,ny_fine,bx,by,s+cs,t+ct);
         mirror_bounds_2d(f1,nx_fine,ny_fine,bx,by);
         view_backward_registration(f0,f1,f1_warp,z,s,t,nx_fine,ny_fine,bx,by,hx_fine,hy_fine);
         lf_4d_setview_with_boundary(lf_res_warp,f1_warp,sx,sy,nx_fine,ny_fine,bx,by,s+cs,t+ct);

      }
   /*don't forget to copy view center.*/
   for( i=0; i<nx_fine;i++)
      for( j=0; j<ny_fine; j++)
         lf_res_warp[cs][ct][i][j] = f0[i+bx][j+by];

   FREE_2D_MATRIX(f0,nx_fine+2*bx,ny_fine+2*by);
   FREE_2D_MATRIX(f1,nx_fine+2*bx,ny_fine+2*by);
   FREE_2D_MATRIX(f1_warp,nx_fine+2*bx,ny_fine+2*by);

}

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

/* prepare memory before calling lf_4d_resample */
lf_4d_warp_interface(
   double* iLF, /* in: 4D LF */
   double* oLF,     /* out: 4D LF */
   double* z,
   int sx,         /* in: origin dim */
   int sy,         /* in: origin dim */
   int nx,
   int ny,
   double hx,        /* out: target dim */
   double hy         /* out: target dim */
)
{

float**** iLF4D;
float**** oLF4D;
float** dis;
float** tmp;
int i,j;
int bx,by;

bx =1;
by = 1;

/* alloc mem */
ALLOC_4D_LF(&iLF4D,sx,sy,nx,ny);
ALLOC_4D_LF(&oLF4D,sx,sy,nx,ny);
ALLOC_2D_MATRIX(&dis,nx+2*bx,ny+2*by);

/*set tmp to zeros*/
/*for(i =0; i<nx; i++)
   for(j =0; j<ny; j++){
      tmp[i][j] = 0.0;
   }
*/
img_2D_from_1D(dis,z,nx,ny,bx,by);

/* copy input */
lf_4d_from_1D(iLF4D,iLF,sx,sy,nx,ny);

/* do resample */
lf_backward_registration(iLF4D,oLF4D,dis,
                        sx,sy,nx,ny,bx,by,hx,hy);
/*copy result back */
lf_4d_to_1D(oLF,oLF4D,sx,sy,nx,ny);

/* free mem */
FREE_4D_LF(iLF4D,sx,sy,nx,ny);
FREE_4D_LF(oLF4D,sx,sy,nx,ny);
FREE_2D_MATRIX(dis,nx+2*bx,ny+2*by);

}


/*lf_backward_registration(lf_res,lf_res_warp,z,
                         sx,sy,nx_fine,ny_fine,bx,by,hx_fine,hy_fine);*/

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


    double* origin_dim;
    double* scales;
    double hx,hy;

    double* iLF;
    double* oLF;
    double* z;


    /* get the value of the scalar input  */
    origin_dim = mxGetPr(piArg[0]);

    sx = (int)origin_dim[0];
    sy = (int)origin_dim[1];
    nx = (int)origin_dim[2];
    ny = (int)origin_dim[3];

    scales = mxGetPr(piArg[1]);

    hx = scales[0];
    hy = scales[1];

    iLF = mxGetPr(piArg[2]);;
    z = mxGetPr(piArg[3]);;

    mexPrintf("Got LF dims: sx %d sy %d nx %d ny %d \n",sx,sy, nx,ny);

    /* get a pointer to the real data in the output matrix */
    /* create the output matrix */
    mwSize dims[4] = {sy, sx, ny,nx};
    poArg[0] =  mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);

    oLF = mxGetPr(poArg[0]);

    lf_4d_warp_interface(iLF,oLF,z,sx,sy,nx,ny,hx,hy);
}
