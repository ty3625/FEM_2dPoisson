/*
 * 2d_femtest.h
 *
 *  Created on: 2016/05/28
 *      Author: TY
 */

#ifndef FEMTEST_H_
#define FEMTEST_H_

int initconstmesh(MESH *msh);
double calcdetA(MESH *msh, int index);
int set_ith_matrix_uxvxuyvy(MESH *msh, int index, double DD, double *M);
int set_ith_vector_f(MESH *msh, int index, double fbar, double *f);
int set_dirichret_bc_bylabel(MESH *msh, int label, double value, double *M, double *f);
int vector_set_zeros(int size, double *u);



#endif /* 2D_FEMTEST_H_ */
