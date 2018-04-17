/*****************************************************************************/
/*
/* 16.10.31 Trung-Hieu Tran
/* PAS - IPVS - Universitaet Stuttgart.
/*
/*****************************************************************************/

#ifndef __COMMONS_H
#define __COMMONS_H

void data_2D_from_1D(float** iIMG2D,double* iIMG,int nx,int ny,int bx,int by);
void data_2D_to_1D(double* oIMG,float** oIMG2D,int nx,int ny,int bx,int by);
void set_bounds_2d(float** out,int nx,int ny,int bx,int by,float val);
void set_val_2d(float** out,int nx,int ny,int bx,int by,double val);

#endif
