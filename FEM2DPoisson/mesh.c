#include <stdio.h>
#include <string.h>
#include "mesh.h"

#define MAXLINE 32768

int findkey(FILE *fp, const char *key)
{
	char temp[MAXLINE];

	do {
	if (NULL == fgets(temp, MAXLINE, fp))
		return 1;
	} while (NULL == strstr(temp, key));

	return 0;
}

int save2dmesh(MESH *msh, const char *filename)
{
	FILE *fp;
	int i;

	fp = fopen(filename, "w");

	fprintf(fp, "Vertices\n");
	fprintf(fp, "%d\n", msh->nvert);
	for (i=0; i<msh->nvert; i++)
		fprintf(fp, "%.5f %.5f %d\n", msh->x[i], msh->y[i], msh->vlabel[i]);
	fprintf(fp, "\n");
	fprintf(fp, "Edegs\n");
	fprintf(fp, "%d\n", msh->nedge);
	for (i=0; i<msh->nedge; i++)
		fprintf(fp, "%d %d %d\n", msh->ep1[i], msh->ep2[i], msh->elabel[i]);
	fprintf(fp, "Triangles\n");
	fprintf(fp, "%d\n", msh->ntriangle);
	for (i=0; i<msh->ntriangle; i++)
		fprintf(fp, "%d %d %d %d\n", msh->tp1[i], msh->tp2[i], msh->tp3[i], msh->tlabel[i]);

	fprintf(fp, "\n");

	fclose (fp);
	return 0;
}

int read2dmesh(MESH *msh, const char *filename)
{
	int i;
	FILE *fp;
	fp = fopen(filename, "r");

	if (NULL == fp)
		return 2;

	if (findkey(fp, "Vertices"))
		goto err;
	fscanf(fp, "%d\n", &(msh->nvert));

	msh->x = (double *)malloc(sizeof(double) * msh->nvert);
	msh->y = (double *)malloc(sizeof(double) * msh->nvert);
	msh->vlabel = (int *)malloc(sizeof(int)  * msh->nvert);

	for (i=0; i<msh->nvert; i++)
		fscanf(fp, "%lf %lf %d\n", msh->x+i, msh->y+i, msh->vlabel+i);

	if (findkey(fp, "Edges"))
		goto err;
	fscanf(fp, "%d\n", &(msh->nedge));

	msh->ep1 = (int *)malloc(sizeof(int) * msh->nedge);
	msh->ep2 = (int *)malloc(sizeof(int) * msh->nedge);
	msh->elabel = (int *)malloc(sizeof(int) * msh->nedge);

	for (i=0; i<msh->nedge; i++)
		fscanf(fp, "%d %d %d\n", msh->ep1+i, msh->ep2+i, msh->elabel+i);

	if (findkey(fp, "Triangles"))
		goto err;
	fscanf(fp, "%d\n", &(msh->ntriangle));

	msh->tp1 = (int *)malloc(sizeof(int) * msh->ntriangle);
	msh->tp2 = (int *)malloc(sizeof(int) * msh->ntriangle);
	msh->tp3 = (int *)malloc(sizeof(int) * msh->ntriangle);
	msh->tlabel = (int *)malloc(sizeof(int) * msh->ntriangle);

	for (i=0; i<msh->ntriangle; i++)
		fscanf(fp, "%d %d %d %d\n", msh->tp1+i, msh->tp2+i, msh->tp3+i, msh->tlabel+i);



	fclose(fp);

	msh->nbvert = 0;
	msh->bv = NULL;

	return 0;

err:
	fclose(fp);
	return 1;
}

