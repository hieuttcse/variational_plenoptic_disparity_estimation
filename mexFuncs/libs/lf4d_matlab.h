/*****************************************************************************/
/*
/* 16.10.31 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#ifndef __LF4D_MATLAB_H
#define __LF4D_MATLAB_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void ALLOC_4D_LF(float***** lf,int sx, int sy, int nx, int ny);
void FREE_4D_LF(float**** lf,int sx, int sy, int nx, int ny);
void lf_4d_extractview(float** view, float**** lf,int sx, int sy, int nx, int ny, int vx, int vy);
void lf_4d_extractview_with_boundary(float** view, float**** lf,int sx, int sy, int nx, int ny, int bx, int by, int vx, int vy);
void lf_4d_updateview(float**** lf,float** view, int sx, int sy, int nx, int ny, int vx, int vy);
void lf_4d_from_5d(float**** lf4d, float***** lf5d, int sx, int sy, int nx, int ny);

void lf_4d_from_1D(float**** lf4d, double * lf1d, int sx, int sy, int nx, int ny);
void lf_4d_to_1D(double * lf1d, float**** lf4d, int sx, int sy, int nx, int ny);


void lf_4d_resample(float**** lf,int sx,int sy,int nx,int ny,float**** lf_res,int nx_fine,int ny_fine,float** tmp);
void lf_4d_setview_with_boundary(float**** lf, float** view, int sx, int sy, int nx, int ny,int bx, int by, int vx, int vy);

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
);

#endif
