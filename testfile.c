#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <rot_trans.h>

int main (int argc, char *argv[])
{
	int i, size=0;
	float x,y,z;
	gsl_matrix *data1;
	gsl_matrix *data2;
	FILE *fp = NULL;
	Transformation T1; 

	if((fp = fopen("./Data/data","r")) == NULL)
		{
			printf("Data file not found.\n");
			return -1;
		}

	fscanf(fp,"%d",&size);
	printf("Number of Data points: %d\n",size);

	data1 = gsl_matrix_alloc (size, DIM);
	data2 = gsl_matrix_alloc (size, DIM);

	for (i = 0; i < 2*size ; ++i)
	{
		fscanf(fp," %f %f %f ", &x, &y, &z);

		if (i < size)
		{
			gsl_matrix_set (data1, i, 0, x);
			gsl_matrix_set (data1, i, 1, y);
			gsl_matrix_set (data1, i, 2, z);
		}
		else 
		{
			gsl_matrix_set (data2, i - size, 0, x);
			gsl_matrix_set (data2, i - size, 1, y);
			gsl_matrix_set (data2, i - size, 2, z);
		}
	}
	/*
		get sample data from Matlab.
	*/

	display(data1,size,DIM, "Data Point 1");
	display(data2,size,DIM, "Data Point 2");

	T1 = findRotationAndTranslation(data1, data2, size);
	
	display(T1.rotationMatrix,DIM,DIM,"Rotation Matrix");
	display(T1.translationVector,DIM,1,"Rotation Matrix");

	gsl_matrix_free(data1);
	gsl_matrix_free(data2);
	return 0;
}