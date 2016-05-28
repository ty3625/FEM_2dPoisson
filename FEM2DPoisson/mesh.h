#ifndef MESH_H_
#define MESH_H_

typedef struct mesh {
	int nvert;
	double *x;
	double *y;
	int *vlabel;
	int ntriangle;
	int *tp1;
	int *tp2;
	int *tp3;
	int *tlabel;
	int nbvert;
	int *bv;
	int nedge;
	int *ep1;
	int *ep2;
	int *elabel;
} MESH;

int read2dmesh(MESH *, const char *);
int save2dmesh(MESH *, const char *);

#endif
