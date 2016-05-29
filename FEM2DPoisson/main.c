// solver for 2D Poisson Problem
// div (grad (u)) + f = 0,   u=u0 at y=0
//
// TY <verne.ty3625@gmail.com>
//
// input: medit .mesh format
//


#include <stdio.h>
#include <stdlib.h>
#include <mkl_lapack.h>
#include "mesh.h"
#include "2d_femtest.h"

#define TGV 1e34

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

	// coefficients
	// A div(grad(u)) + f0 = 0, u=u0 at y==0
	const double A = 1.0;
	const double f0 = 1.0;
	const double u0 = 2.3;


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

	// generate Stiffness matrix
	// A * int2d( ux vx + uy vy ) = int2d (f0 * v)
	for (i=0; i<msh->ntriangle; i++) {
		set_ith_matrix_uxvxuyvy(msh, i, A, M);
		set_ith_vector_f(msh, i, f0, f);
	}

//	set_dirichret_bc_bylabel(msh, 1, 2.3, M, f);  // dirichret boundary condition

// dirichret boundary condition for meshes y == 0.0
// u0 = 2.3;
	for (i=0; i<msh->nvert; i++) {
		if (msh->y[i] < 0.01) {
			M[i + i * nvert] += TGV;
			f[i] += TGV * u0;
		}
	}

	printf("Stiffness Matrix M\n");
	for (j=0; j<nvert; j++) {
		for (i=0; i<nvert; i++) {
			printf("%.3e\t", M[j + i * nvert]);
		}
		printf("\n");
	}
	printf("\n");

	printf("Load Vector\n");
	for (i=0; i<nvert; i++)
		printf("%.3e\n", f[i]);
	printf("\n");
	printf("\n");


	printf("solve...\n");
	dsytrf(uplo, &n, M, &lda, p, work, &lwork, &info);
	printf("dsytrf: %d\n", info);

	dsytrs(uplo, &n, &nrhs, M, &lda, p, f, &ldb, &info);
	printf("dsytrs: %d\n", info);
	printf("\n");

	printf("Displacements\n");
	for (i=0; i<nvert; i++)
		printf("%.12e\n", f[i]);
	printf("\n");


	return 0;
}
