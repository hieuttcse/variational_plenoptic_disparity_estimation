/*****************************************************************************/
/*
/* 16.10.27 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#include "mex.h"
#include "math.h"



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





/* ------------------------------------------------------------------------- */

void horn_schunck_warp_sor_lf
(
        /*****************************************************/
        float **P,   /* in     : entry 11 of the motion tensor            */
        float **Q,   /* in     : entry 33 of the motion tensor            */
        float **W,
        float **dz,     /* in+out : x-component of flow increment            */
        float **z,      /* in     : x-component of flow field                */
        int nx,         /* in     : size in x-direction                      */
        int ny,         /* in     : size in y-direction                      */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx,       /* in     : grid spacing in x-direction              */
        float hy,       /* in     : grid spacing in y-direction              */
        float alpha,    /* in     : smoothness weight                        */
        float sor     /* in     : SOR overrelaxation parameter             */
                        /*****************************************************/
)

/*
   Computes one SOR iteration
 */

{
        /*****************************************************/
        int i,j;        /* loop variables                                    */
        float hx_2,hy_2; /* time saver variables                              */
        float xp,xm,yp,ym; /* neighbourhood weights                             */
        float sum;      /* central weight                                    */
        float p,q;
        float wxm, wxp, wym,wyp;
        float wx;
                        /*****************************************************/

/* define time saver variables */
        hx_2=alpha/(hx*hx);
        hy_2=alpha/(hy*hy);


/* set boundaries zero */
        /*set_bounds_2d(z,nx,ny,bx,by,0.0);
        set_bounds_2d(dz,nx,ny,bx,by,0.0);
*/

        for(i=bx; i<nx+bx; i++)
                for(j=by; j<ny+by; j++)
                {
                        /* compute weights */

                        wxm = (W[i-1][j] + W[i][j])/2;
                        wxp = (W[i+1][j] + W[i][j])/2;
                        wym = (W[i][j-1] + W[i][j])/2;
                        wyp = (W[i][j+1] + W[i][j])/2;

                        xp =  (i<nx+bx-1) * hx_2*wxp;
                        xm =  (i>bx)      * hx_2*wxm;
                        yp =  (j<ny+by-1) * hy_2*wyp;
                        ym =  (j>by)      * hy_2*wym;




                        p = P[i][j];
                        q = Q[i][j];

                        sum = xp + xm + yp + ym;

                        /* perform SOR iteration */

                        dz[i][j] = (1-sor) * dz[i][j] + sor *
                                   ( ( q
                                       - xm*( dz[i-1][j] + z[i-1][j] )
                                       - xp*( dz[i+1][j] + z[i+1][j] )
                                       - ym*( dz[i][j-1] + z[i][j-1] )
                                       - yp*( dz[i][j+1] + z[i][j+1] )
                                       + sum * ( z[i][j] ))
                                     /(-p -sum ) );

                        if(dz[i][j] != dz[i][j]){
                           z[i][j] = 1;
                        }

                        /* --------------------------------------- */

                }
}

void update_robust_weight
(
   float** P, /* remember P, Q, dz size nx+2 x ny+2 */
   float** Q,
   double* Jg_11,
   double* Jg_22,
   double* Jg_12,
   double* JG_11,
   double* JG_22,
   double* JG_12,
   int nx,
   int ny,
   double gamma, /* weight for gradient constancy */
   double epsi,
   float** dz
){
   int i,j,s,t;
   float wg, wG;
   float dis;
   float jg11,jg12,jg22;
   float jG11,jG12,jG22;

   for( i =0; i<nx; i++){
      for(j=0; j<ny; j++){
         dis = dz[i+1][j+1];
         jg11 = Jg_11[i*ny+j];
         jg12 = Jg_12[i*ny+j];
         jg22 = Jg_22[i*ny+j];
         wg = jg11*dis*dis + 2*dis*jg12 + jg22;
         if(wg < 0)
            wg = 0;
         wg = 1.0 / sqrt(wg + epsi);

         jG11 = JG_11[i*ny+j];
         jG12 = JG_12[i*ny+j];
         jG22 = JG_22[i*ny+j];
         wG = jG11*dis*dis + 2*dis*jG12 + jG22;
         if(wG <0)
            wG =0;
         wG = 1.0 / sqrt(wG + epsi);


         P[i+1][j+1] = wg*jg11 + wG*jG11*gamma;
         Q[i+1][j+1] = wg*jg12 + wG*jG12*gamma;
/*
         if(i ==147 && j==142){
            P[i+1][j+1]=0;
         }
*/
      /*   if(P[i+1][j+1] != P[i+1][j+1] || Q[i+1][j+1]!=Q[i+1][j+1]){
            P[i+1][j+1]=0;
         }*/
      }
   }
}



void update_flow_driven_weight
(
   float** W,
   int nx,
   int ny,
   double epsi,
   float** dz,
   float** z
)
{
   float** zx;
   float** zy;
   int i;
   int j;

   float** dzx;
   float** dzy;

   float w;


   ALLOC_2D_MATRIX(&zx,nx+2,ny+2);
   ALLOC_2D_MATRIX(&zy,nx+2,ny+2);
   ALLOC_2D_MATRIX(&dzx,nx+2,ny+2);
   ALLOC_2D_MATRIX(&dzy,nx+2,ny+2);


   /* compute first derivatev of z, dz */
   for( i =1; i<nx+1; i++){
      for( j=1; j<ny+1; j++){
         zx[i][j] = (z[i+1][j] - z[i-1][j])/2.0;
         zy[i][j] = (z[i][j+1] - z[i][j-1])/2.0;
         dzx[i][j] = (dz[i+1][j] - dz[i-1][j])/2.0;
         dzy[i][j] = (dz[i][j+1] - dz[i][j-1])/2.0;

         w = (zx[i][j]+dzx[i][j])*(zx[i][j]+dzx[i][j]) + (zy[i][j] + dzy[i][j])* (zy[i][j]+dzy[i][j]);
         W[i][j] = 1.0 /(2.0*sqrt(w+epsi));
      }
   }

   FREE_2D_MATRIX(zx,nx+2,ny+2);
   FREE_2D_MATRIX(zy,nx+2,ny+2);
   FREE_2D_MATRIX(dzx,nx+2,ny+2);
   FREE_2D_MATRIX(dzy,nx+2,ny+2);
}


/* The gateway function */
/* dz = warp_solver_gs([nx ny bx by],[hx hy],[100 1.5 100],J11,J22,J12,z); */
void mexFunction( int noArg, mxArray *poArg[],
                  int niArg, const mxArray *piArg[])
{
    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    int i;
    int j;

    int sx;
    int sy;
    int nx;
    int ny;

    int iteration;


    float hx;
    float hy;

    double* Jg_11;
    double* Jg_12;
    double* Jg_22;

    double* JG_11;
    double* JG_12;
    double* JG_22;

    double* Zinit;

    double alpha;
    double gamma;
    double epsi;

    double sor;

    double* dimension;
    double* scale;
    double* params;
    double* LF;

    float** z;
    float** dz;

    float** P;
    float** Q;
    float** W;

    /* get the value of the scalar input  */
    dimension = mxGetPr(piArg[0]);
    nx = (int)dimension[0];
    ny = (int)dimension[1];

    scale = mxGetPr(piArg[1]);
    hx = (int)scale[0];
    hy = (int)scale[1];

    params = mxGetPr(piArg[2]);
    alpha = params[0];
    gamma = params[1];
    epsi  = params[2];
    sor = params[3];
    iteration = params[4];

    Jg_11 = mxGetPr(piArg[3]);
    Jg_22 = mxGetPr(piArg[4]);
    Jg_12 = mxGetPr(piArg[5]);

    JG_11 = mxGetPr(piArg[6]);
    JG_22 = mxGetPr(piArg[7]);
    JG_12 = mxGetPr(piArg[8]);

    Zinit = mxGetPr(piArg[9]);


    mexPrintf("Got dim %d %d \n",nx,ny);
    mexPrintf("Got scaling factor %f %f \n",hx,hy);

    mexPrintf("Got params alpha %f gamma %f epsi %f sor %f \n",alpha,gamma, epsi, sor);

    /*    mexPrintf("Got J_11, 1 %f 2 %f 5 %f\n",J_11[1], J_11[2], J_11[5]);*/

   /* ALLOC MEM */
   ALLOC_2D_MATRIX(&z,nx+2,ny+2);
   ALLOC_2D_MATRIX(&dz,nx+2,ny+2);
   ALLOC_2D_MATRIX(&P,nx+2,ny+2);
   ALLOC_2D_MATRIX(&Q,nx+2,ny+2);
      ALLOC_2D_MATRIX(&W,nx+2,ny+2);

    /* get a pointer to the real data in the output matrix */
    /* create the output matrix */
    poArg[0] =  mxCreateDoubleMatrix(ny,nx,mxREAL);
    outMatrix = mxGetPr(poArg[0]);

    /* create compute filed. */
    /* copy input to 2D array */
    img_2D_from_1D(z,Zinit,nx,ny,1,1);

    set_val_2d(dz,nx,ny,1,1,0.0);
    set_val_2d(P,nx,ny,1,1,0.0);
    set_val_2d(Q,nx,ny,1,1,0.0);
    set_val_2d(W,nx,ny,1,1,0.0);
    set_bounds_2d(z,nx,ny,1,1,0.0);


    for(i = 0; i<iteration; i++){
      set_bounds_2d(dz,nx,ny,1,1,0.0);
      set_bounds_2d(W,nx,ny,1,1,0.0);
      update_robust_weight(P,Q,Jg_11,Jg_22,Jg_12, JG_11,JG_22,JG_12, nx,ny, gamma, epsi, dz);
      update_flow_driven_weight(W,nx,ny, epsi, dz, z);
      horn_schunck_warp_sor_lf(P,Q,W, dz,z,nx,ny,1,1,hx,hy,alpha,sor);
    }

    /* copy result back */
    img_2D_to_1D(outMatrix,dz,nx,ny,1,1);

    /* FREE MEM */
    FREE_2D_MATRIX(z,nx+2,ny+2);
    FREE_2D_MATRIX(dz,nx+2,ny+2);
    FREE_2D_MATRIX(P,nx+2,ny+2);
    FREE_2D_MATRIX(Q,nx+2,ny+2);
    FREE_2D_MATRIX(W,nx+2,ny+2);
}
