/*****************************************************************************/
/*
/* 29.10.31 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

/*** Basic functions should be covered ***/
/*
1. create 4D LF from 1D colored-LF
2. extract abitrary view.
*/

#include"lf4d_matlab.h"


float sRGB_to_linear(float x){
	if(x<0.04045) return x/12.92;
	return pow((x+0.055)/1.055,2.4);
}
/* Apply gramma correction*/
float linear_to_sRGB(float y){
	if(y<=0.0031308) return 19.92*y;
	return 1.055*pow(y,1/2.4)-0.055;
}
float convertRGBtoGray(float* colors){
	float R,G,B,est;
	R = sRGB_to_linear(colors[0]);
	G = sRGB_to_linear(colors[1]);
	B = sRGB_to_linear(colors[2]);
	est = 0.2126 * R + 0.7152 * G + 0.0722*B;
	return linear_to_sRGB(est);
}



void ALLOC_4D_LF(float***** lf,int sx, int sy, int nx, int ny){
   float**** array;
	int s,t,i,j;
    array = (float****)malloc( sizeof(float***) * sx );
   for( s =0; s<sx; s++){
      array[s] = (float***)malloc( sizeof(float**) * sy);
      for( t=0; t<sy; t++){
         array[s][t] = (float**)malloc( sizeof(float*) * nx);
         for( i=0; i<nx; i++)
            array[s][t][i] = (float*)malloc( sizeof(float) * ny);
      }
   }
   *lf = array;
}

void FREE_4D_LF(float**** lf,int sx, int sy, int nx, int ny){
	int s,t,i;
   for( s=0; s<sx; s++){
      for( t=0; t<sy; t++){
         for( i=0; i<sx; i++){
               free(lf[s][t][i]);
         }
         free(lf[s][t]);
      }
      free(lf[s]);
   }
   free(lf);
}



/* ------------------------------------------------------------------------- */


void my_resample_1d

(
                /*************************************************************/
    float *u,   /* in     : input vector, size 1..n                          */
    int   n,    /* in     : size of input vector                             */
    int   m,    /* in     : size of output vector                            */
    float *v    /* out    : output vector, size 1..m                         */
                /*************************************************************/
    )

/* Area-based resampling: Transforms a 1D image u of size n into an image v  */
/* of size m by integration over piecewise constant functions. Conservative. */

{
                         /****************************************************/
int     i, k;            /* loop variables                                   */
float   hu, hv;          /* grid sizes                                       */
float   uleft, uright;   /* boundaries                                       */
float   vleft, vright;   /* boundaries                                       */
float   fac;             /* normalization factor                             */
                         /****************************************************/

/*****************************************************************************/
/* (1/2) Special cases of area-based resampling are computed efficiently     */
/*****************************************************************************/

/* fast interpolation for output images of even size */
if (m==2*n)
 {
     /* one cell is devided in two cells with equal value */
     for (i=1; i<=n; i++)
     {
     v[i*2-1]=u[i];
     v[i*2  ]=u[i];
     }
     return;
 }


/* fast restriction for input images of even size */
if (2*m==n)
 {
     /* two celss are melted to a larger cell with averaged value */
     for (i=1; i<=m; i++)
     {
     v[i]=0.5*(u[i*2-1]+u[i*2]);
     }

     return;
     }


/*****************************************************************************/
/* (2/2) Remaining cases require more complex algorithm                      */
/*****************************************************************************/


/* initializations */
                            /*************************************************/
hu    = 1.0 / (float)n;     /* grid size of u                                */
hv    = 1.0 / (float)m;     /* grid size of v                                */
uleft = 0.0;                /* left interval boundary of u                   */
vleft = 0.0;                /* left interval boundary of v                   */
k     = 1;                  /* index for u                                   */
fac   = hu / hv;            /* for normalization                             */
                            /*************************************************/


/*---- loop ----*/

for (i=1; i<=m; i++)
    /* calculate v[i] by integrating the piecewise constant function u */
    {
    /* calculate right interval boundaries */
    uright = uleft + hu;
    vright = vleft + hv;

    if (uright > vright)
       /* since uleft <= vleft, the entire v-cell i is in the u-cell k */
       v[i] = u[k];
    else
       {
       /* consider fraction alpha of the u-cell k in v-cell i */
       v[i] = (uright - vleft) * n * u[k++];

       /* update */
       uright = uright + hu;

       /* consider entire u-cells inside v-cell i */
       while (uright <= vright)
             /* u-cell k lies entirely in v-cell i; sum up */
             {
             v[i] = v[i] + u[k++];
             uright = uright + hu;
             }

       /* consider fraction beta of the u-cell k in v-cell i */
       v[i] = v[i] + (1.0 - (uright - vright) * n) * u[k];

       /* normalization */
       v[i] = v[i] * fac;
       } /* else */

    /* update */
    uleft = uright - hu;
    vleft = vright;
    /* now it holds: uleft <= vleft */
    }  /* for i */


return;

}

void MY_ALLOC_VECTOR(float** out, int size){
	float* array = (float*) malloc(sizeof(float)*size);
	*out = array;
}

void MY_FREE_VECTOR(float* out, int size){
	free(out);
}

void my_resample_2d_x

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    mx       /* in   : x dimension of output image                        */
                /*************************************************************/
)

/* resample a 2-D image in x-direction using area-based resampling */

{
                         /****************************************************/
int    i, j;             /* loop variables                                   */
float  *uhelp, *vhelp;   /* auxiliary vectors                                */
                         /****************************************************/


/* allocate memory */
MY_ALLOC_VECTOR ( &uhelp, nx+2);
MY_ALLOC_VECTOR ( &vhelp, mx+2);

/* resample image linewise in x-direction */
for (j=by; j<ny+by; j++)
 {
     /* initialise left boundary of 1-D array with zero */
     uhelp[0]=u[bx][j];

     /* copy current line in this 1-D array */
     for (i=bx; i<nx+bx; i++)
	 uhelp[i-bx+1] = u[i][j];

    /* initialise right boundary of this 1-D array with zero */
    uhelp[nx+1]=u[nx+bx-1][j];

    /* resample this 1-D array */
    my_resample_1d (uhelp, nx, mx, vhelp);

    /* copy resmapled array in corresponding output line */
    for (i=bx; i<mx+bx; i++)
		u_out[i][j] = vhelp[i-bx+1];
    }

/* free memory */
MY_FREE_VECTOR (uhelp, nx+2);
MY_FREE_VECTOR (vhelp, mx+2);

return;
}

void my_resample_2d_y

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    my       /* in   : y dimension of output image                        */
                /*************************************************************/
)

/* resample a 2-D image in y-direction using area-based resampling */

{
                         /****************************************************/
int    i, j;             /* loop variables                                   */
float  *uhelp, *vhelp;   /* auxiliary vectors                                */
                         /****************************************************/


/* allocate memory */
MY_ALLOC_VECTOR (&uhelp, ny+2);
MY_ALLOC_VECTOR (&vhelp, my+2);

/* resample image columnwise in y-direction */
for (i=bx; i<nx+bx; i++)
 {
     /* initialsie left boundary of 1-D array with zero */
     uhelp[0]=u[i][by];

     /* copy current column in this 1-D array */
     for (j=by; j<ny+by; j++)
	 uhelp[j-by+1] = u[i][j];

     /* initialise right boundary of this 1-D array with zero */
     uhelp[ny+1]=u[i][ny+by-1];

     /* resample this 1-D array */
     my_resample_1d (uhelp, ny, my, vhelp);

     /* copy resmapled array in corresponding output column */
     for (j=by; j<my+by; j++)
	 u_out[i][j] = vhelp[j-by+1];

 }

/* free memory */
MY_FREE_VECTOR (uhelp, ny+2);
MY_FREE_VECTOR (vhelp, my+2);

return;
}


/*--------------------------------------------------------------------------*/



void my_resample_2d

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    mx,      /* in   : x dimension of output image                        */
int    my,      /* in   : y dimension of output image                        */
float  **tmp    /* tmp  : temporary array of size (nx * my)                  */
                /*************************************************************/
)

/* resample a 2-D image using area-based resampling */

{

    /* if interpolation */
    if (my>=ny)
    {
	/* resample first in x-direction */
	my_resample_2d_x(u, nx, ny, bx, by, tmp, mx);
	/* resample then in y-direction */
	my_resample_2d_y(tmp, mx, ny, bx, by, u_out, my);
    }
    /* if restriction */
    else
    {
	/*  resample first in y-direction */
	my_resample_2d_y(u, nx, ny, bx, by, tmp, my);
        /*  resample tehn in x-direction */
	my_resample_2d_x(tmp, nx, my, bx, by, u_out, mx);
    }

return;

}

void lf_4d_resample(float**** lf,int sx,int sy,int nx,int ny,float**** lf_res,int nx_fine,int ny_fine,float** tmp){
	int s,t;
	for( s=0; s<sx; s++)
		for( t=0; t<sy; t++){
			my_resample_2d(lf[s][t], nx,ny,0,0,lf_res[s][t],nx_fine,ny_fine,tmp);
		}
}

void lf_4d_extractview(float** view, float**** lf, int sx, int sy, int nx, int ny, int vx, int vy){
	int i,j;

   if(vx <0 || vx>=sx||vy<0 || vy>=sy){
      fprintf(stderr," EXTRACT VIEW: out of view index\n");
      fprintf(stderr," EXTRACT VIEW: extract %d x %d out of %d x %d x %d x %d \n", vx,vy, sx,sy,nx,ny);
      exit(0);
   }
   for( i =0; i<nx; i++)
      for( j=0; j<ny; j++)
         view[i][j] = lf[vx][vy][i][j];
}

void lf_4d_extractview_with_boundary(float** view, float**** lf, int sx, int sy, int nx, int ny,int bx, int by, int vx, int vy){
	int i;
	int j;

   if(vx <0 || vx>=sx||vy<0 || vy>=sy){
      fprintf(stderr," EXTRACT VIEW: out of view index\n");
      fprintf(stderr," EXTRACT VIEW: extract %d x %d out of %d x %d x %d x %d \n", vx,vy, sx,sy,nx,ny);
      exit(0);
   }
   for( i =0; i<nx; i++)
      for( j=0; j<ny; j++)
         view[i+bx][j+by] = lf[vx][vy][i][j];
}

void lf_4d_setview_with_boundary(float**** lf, float** view, int sx, int sy, int nx, int ny,int bx, int by, int vx, int vy){
	int i,j;

   if(vx <0 || vx>=sx||vy<0 || vy>=sy){
      fprintf(stderr," SET VIEW: out of view index\n");
      fprintf(stderr," SET VIEW: extract %d x %d out of %d x %d x %d x %d \n", vx,vy, sx,sy,nx,ny);
      exit(0);
   }
   for( i =0; i<nx; i++)
      for( j=0; j<ny; j++)
      	lf[vx][vy][i][j] = view[i+bx][j+by];
}


void lf_4d_updateview(float**** lf,float** view, int sx, int sy, int nx, int ny, int vx, int vy){
	int i,j;

   if(vx <0 || vx>=sx||vy<0 || vy>=sy){
      fprintf(stderr," UPDATE VIEW: out of view index\n");
      fprintf(stderr," UPDATE VIEW: extract %d x %d out of %d x %d x %d x %d \n", vx,vy, sx,sy,nx,ny);
      exit(0);
   }
   for( i =0; i<nx; i++)
      for( j=0; j<ny; j++)
         lf[vx][vy][i][j]= view[i][j];
}

/**
* Convert from 5d lf to 4d lf
**/
void lf_4d_from_5d(float**** lf4d, float***** lf5d, int sx, int sy, int nx, int ny){
	int s,t,i,j;

   for( s=0; s<sx; s++){
      for( t=0; t<sy; t++){
         for( i=0; i<nx; i++){
            for( j=0; j<ny; j++){
               lf4d[s][t][i][j] = convertRGBtoGray(lf5d[s][t][i][j]);
            }
         }
      }
   }
}


void lf_4d_from_1D(float**** lf4d, double * lf1d, int sx, int sy, int nx, int ny){
	int s,t,i,j;

   for( s=0; s<sx; s++)
      for( t=0; t<sy; t++)
         for( i =0; i<nx; i++)
            for( j=0; j<ny; j++)
                  lf4d[s][t][i][j]= (float)lf1d[i*ny*sx*sy + j*sx*sy +s*sy + t];
}


void lf_4d_to_1D(double* lf1d, float**** lf4d, int sx, int sy, int nx, int ny){
	int s,t,i,j;

   for( s=0; s<sx; s++)
      for( t=0; t<sy; t++)
         for( i =0; i<nx; i++)
            for( j=0; j<ny; j++)
                  lf1d[i*ny*sx*sy + j*sx*sy +s*sy + t] = lf4d[s][t][i][j];
}
