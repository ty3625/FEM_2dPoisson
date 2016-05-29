#include <stdio.h>
#include <mkl_blas.h>
#include "mesh.h"
#include "2d_femtest.h"

int initconstmesh(MESH *msh)
{
	int nvert = 9;
	static double x[] = {0, 1, 2, 0 ,1, 2, 0, 1, 2};
	static double y[] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
	static int vlabel[] = {1, 1, 1, 0, 0, 0, 3, 3, 3};

	int ntriangle = 8;
	static int tp1[] = {1, 1, 4, 7, 2, 3, 5, 9};
	static int tp2[] = {2, 5, 5, 5, 3, 6, 6, 8};
	static int tp3[] = {5, 4, 7, 8, 5, 5, 9, 5};

	int nbvert = 8;
	static int bv[] = {1, 2, 3, 6, 9, 8, 7, 4};

	msh->nvert = nvert;
	msh->x = x;
	msh->y = y;
	msh->vlabel = vlabel;
	msh->ntriangle = ntriangle;
	msh->tp1 = tp1;
	msh->tp2 = tp2;
	msh->tp3 = tp3;
	msh->tlabel = NULL;
	msh->nbvert = nbvert;
	msh->bv = bv;

	return 0;
}

double calcdetA(MESH *msh, int index) {
	double *x, *y;
	int *tp1, *tp2, *tp3;
	x = msh->x;
	y = msh->y;
	tp1 = msh->tp1;
	tp2 = msh->tp2;
	tp3 = msh->tp3;

	return
	  x[tp2[index]-1] * y[tp3[index]-1]
	+ x[tp1[index]-1] * y[tp2[index]-1]
	+ x[tp3[index]-1] * y[tp1[index]-1]
	- x[tp1[index]-1] * y[tp3[index]-1]
	- x[tp2[index]-1] * y[tp1[index]-1]
	- x[tp3[index]-1] * y[tp2[index]-1];
}

int set_ith_matrix_uxvxuyvy(MESH *msh, int index, double DD, double *M) {
	double Be[6];
	double Be2[9];
	double detA;
	int lda = 2;
	int ldc = 3;
	double alpha, beta;
	int nrow;
	int ncol;

	double *x, *y;
	int *tp1, *tp2, *tp3;

	x = msh->x;
	y = msh->y;
	tp1 = msh->tp1;
	tp2 = msh->tp2;
	tp3 = msh->tp3;

	nrow = ncol = msh->nvert;

	// int2d( ux vx + uy vy )

	detA = calcdetA(msh, index);

	Be[0] = y[tp2[index]-1] - y[tp3[index]-1];
	Be[2] = y[tp3[index]-1] - y[tp1[index]-1];
	Be[4] = y[tp1[index]-1] - y[tp2[index]-1];
	Be[1] = x[tp3[index]-1] - x[tp2[index]-1];
	Be[3] = x[tp1[index]-1] - x[tp3[index]-1];
	Be[5] = x[tp2[index]-1] - x[tp1[index]-1];

	alpha = DD / detA / 2.0;
	beta = 0.0;
	dgemm("T", "N", &ldc, &ldc, &lda, &alpha,
		Be, &lda, Be, &lda, &beta, Be2, &ldc);

	M[(tp1[index]-1) + (tp1[index]-1) * nrow] += Be2[0];
	M[(tp2[index]-1) + (tp1[index]-1) * nrow] += Be2[1];
	M[(tp3[index]-1) + (tp1[index]-1) * nrow] += Be2[2];
	M[(tp1[index]-1) + (tp2[index]-1) * nrow] += Be2[3];
	M[(tp2[index]-1) + (tp2[index]-1) * nrow] += Be2[4];
	M[(tp3[index]-1) + (tp2[index]-1) * nrow] += Be2[5];
	M[(tp1[index]-1) + (tp3[index]-1) * nrow] += Be2[6];
	M[(tp2[index]-1) + (tp3[index]-1) * nrow] += Be2[7];
	M[(tp3[index]-1) + (tp3[index]-1) * nrow] += Be2[8];

	return 0;
}

int set_ith_vector_f(MESH *msh, int index, double fbar, double *f) {
	double detA;
//	double *x, *y;
	int *tp1, *tp2, *tp3;

//	x = msh->x;
//	y = msh->y;
	tp1 = msh->tp1;
	tp2 = msh->tp2;
	tp3 = msh->tp3;

	// int2d ( fbar v )

	detA = calcdetA(msh, index);

	f[tp1[index]-1] += detA / 6.0 * fbar;
	f[tp2[index]-1] += detA / 6.0 * fbar;
	f[tp3[index]-1] += detA / 6.0 * fbar;

	return 0;
}

int set_dirichret_bc_bylabel(MESH *msh, int label, double value, double *M, double *f)
{
	int i;
	const double tbn = 1e34;

	for (i=0; i<msh->nvert; i++) {	// dirichret boundary condition
		if (msh->vlabel[i] != label)
			continue;
		M[i + msh->nvert * i] += tbn;
		f[i] += tbn * value;
	}
	return 0;
}

int vector_set_zeros(int size, double *u)
{
	while (size-->0)
		u[size] = 0.0;
	return 0;
}
