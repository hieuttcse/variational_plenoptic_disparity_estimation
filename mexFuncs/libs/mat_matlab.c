/*****************************************************************************/
/*
/* 16.10.31 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/
#include"mat_matlab.h"

void ALLOC_2D_MATRIX(float*** mat, int nx, int ny){
   int i;
   float** array;
   array = (float**)malloc( sizeof(float*) * nx);
   for( i=0; i<nx; i++)
      array[i] = (float*)malloc( sizeof(float) * ny);
   *mat = array;
}

void FREE_2D_MATRIX(float** mat,int nx, int ny){
   int i;

   for( i=0; i<nx; i++){
         free(mat[i]);
   }
   free(mat);
}


void ALLOC_3D_MATRIX(float**** mat, int nx, int ny, int nc){
   int i,j;
   float*** array;

   array = (float***)malloc( sizeof(float**) * nx);
   for( i=0; i<nx; i++){
      array[i] = (float**)malloc( sizeof(float*) * ny);
      for( j=0; j<ny; j++){
         array[i][j] = (float*)malloc( sizeof(float) * nc);
      }
   }
   *mat = array;
}

void FREE_3D_MATRIX(float*** mat,int nx, int ny, int nc){
   int i,j;

   for( i=0; i<nx; i++){
      for( j=0; j<ny;j++)
         free(mat[i][j]);
      free(mat[i]);
   }
   free(mat);
}
