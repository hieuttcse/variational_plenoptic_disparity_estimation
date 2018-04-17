/*****************************************************************************/
/*
/* 16.10.31 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/


#ifndef __MAT_MATLAB_H
#define __MAT_MATLAB_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void ALLOC_2D_MATRIX(float*** mat, int nx, int ny);

void FREE_2D_MATRIX(float** mat,int nx, int ny);

void ALLOC_3D_MATRIX(float**** mat, int nx, int ny, int nc);

void FREE_3D_MATRIX(float*** mat,int nx, int ny, int nc);

#endif
