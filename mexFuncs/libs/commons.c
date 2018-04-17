#include"commons.h"



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



void data_2D_from_1D(float** iIMG2D,double* iIMG,int nx,int ny,int bx,int by){
   int i,j;
   int H = ny+2*by;

   for(i = bx; i<nx+bx; i++)
      for(j =by; j<ny+by; j++){
         iIMG2D[i][j] = (float)iIMG[(i-bx)*ny + j-by];
      }
}

void data_2D_to_1D(double* oIMG,float** oIMG2D,int nx,int ny,int bx,int by){
   int i,j;
   int H = ny+2*by;

   for(i = bx; i<nx+bx; i++)
      for(j =by; j<ny+by; j++){
         oIMG[(i-bx)*ny +(j-by)] = oIMG2D[i][j];
      }
}


void set_bounds_2d(float** out,int nx,int ny,int bx,int by,float val){
   int i;
   int j;
   int W = nx+2*bx;
   int H = ny+2*by;

   for( i =0; i<bx; i++)
      for(j=0; j<H; j++)
         out[i][j] = val;
   for( i =nx+bx; i<W; i++)
      for(j=0; j<H; j++)
         out[i][j] = val;
   for( i =0; i<W; i++)
      for(j=0; j<by; j++)
         out[i][j] = val;
   for( i =0; i<W; i++)
      for(j=ny+by; j<H; j++)
         out[i][j] = val;
}

void set_val_2d(float** out,int nx,int ny,int bx,int by,double val){
   int i;
   int j;
   int W = nx+2*bx;
   int H = ny+2*by;

   for( i =0; i<W; i++)
      for(j=0; j<H; j++)
         out[i][j] = val;
}
