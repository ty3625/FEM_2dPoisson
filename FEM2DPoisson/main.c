#include <stdio.h>
#include <stdlib.h>
#include <mkl_lapack.h>
#include "mesh.h"
#include "2d_femtest.h"

int main(int argc, char **argv)
{
	double *M;
	double *f;
	int i, j;
	int n, lda, lwork, info;
	double *work;
	int *p;
	int nrhs, ldb;
	const char *uplo = "U";
	MESH *msh;
	MESH mesh1;
	int nvert;


	msh = &mesh1;
	initconstmesh(msh);
	read2dmesh(msh, "smallmesh.mesh");
	save2dmesh(msh, "testmesh.txt");

	nvert = msh->nvert;

	M = (double *)malloc(sizeof(double) * 1 * nvert * nvert);
	f = (double *)malloc(sizeof(double) * 1 * nvert);

	n = lda = ldb = lwork = nvert;
	nrhs = 1;
	p = (int *) malloc(sizeof(int) * nvert);
	work = (double *)malloc(sizeof(double) * 1 * lwork);

	vector_set_zeros(1 * nvert * nvert, M);
	vector_set_zeros(1 * nvert, f);

	for (i=0; i<msh->ntriangle; i++) {		// generate Stiffness matrix
		// 1.0 * int2d( ux vx + uy vy ) = int2d (1.0 * v)
		set_ith_matrix_uxvxuyvy(msh, i, 1.0, M);
		set_ith_vector_f(msh, i, 1.0, f);
	}

//	set_dirichret_bc_bylabel(msh, 1, 2.3, M, f);  // dirichret boundary condition
	for (i=0; i<msh->nvert; i++) {
		if (msh->y[i] < 0.01) {
			M[i + i * nvert] += 1e34;
			f[i] += 1e34 * 2.3;
		}
	}

	for (j=0; j<nvert; j++) {
		for (i=0; i<nvert; i++) {
			printf("%.3e\t", M[j + i * nvert]);
		}
		printf("\n");
	}
	for (i=0; i<nvert; i++)
		printf("%.3e\n", f[i]);
	printf("\n");


	printf("solve\n");
	dsytrf(uplo, &n, M, &lda, p, work, &lwork, &info);
	printf("dsytrf: %d\n", info);

	dsytrs(uplo, &n, &nrhs, M, &lda, p, f, &ldb, &info);
	printf("dsytrs: %d\n", info);
	printf("\n");

	for (i=0; i<nvert; i++)
		printf("%.12e\n", f[i]);
	printf("\n");


	return 0;
}
