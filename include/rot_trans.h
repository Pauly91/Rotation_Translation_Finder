#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define DIM 3 // Working in the 3 dimensional space


typedef struct Dim3Points
{
	float x;
	float y;
	float z;
} Dim3Point;


typedef struct Centroids
{
	Dim3Point centroidOfData1;
	Dim3Point centroidOfData2;
} Centroid;


typedef struct Transformations
{
	gsl_matrix *rotationMatrix;
	gsl_matrix *translationVector;  // check if data type vector will workout.
} Transformation;



void display(gsl_matrix *m,int row, int column, char *name)
{
	int i,j;
	float holder;
	printf("\n\n");
	printf(" %s = \n\n",name);
	for (i = 0; i < row; ++i)
	{
		for (j = 0; j < column; ++j)
		{
			holder = gsl_matrix_get (m, i, j);
			printf(" %f ",holder);
		}
		printf("\n");
	}
	printf("\n\n");
}


Centroid findCentroid(gsl_matrix *data1, gsl_matrix *data2, int size)
{
	
	int i;
	float sum_d1_x = 0.0, sum_d1_y = 0.0, sum_d1_z = 0.0;
	float sum_d2_x = 0.0, sum_d2_y = 0.0, sum_d2_z = 0.0;
	float holder;
	Centroid CC;
	for (i = 0; i < size; ++i)
	{
		holder =  gsl_matrix_get (data1, i, 0);
		sum_d1_x += holder;
		holder =  gsl_matrix_get (data1, i, 1);
		sum_d1_y += holder;
		holder =  gsl_matrix_get (data1, i, 2);
		sum_d1_z += holder;

		holder =  gsl_matrix_get (data2, i, 0);
		sum_d2_x += holder;
		holder =  gsl_matrix_get (data2, i, 1);
		sum_d2_y += holder;
		holder =  gsl_matrix_get (data2, i, 2);
		sum_d2_z += holder;

	}

	CC.centroidOfData1.x = sum_d1_x / size;
	CC.centroidOfData1.y = sum_d1_y / size;
	CC.centroidOfData1.z = sum_d1_z / size;

	CC.centroidOfData2.x = sum_d2_x / size;
	CC.centroidOfData2.y = sum_d2_y / size;
	CC.centroidOfData2.z = sum_d2_z / size;


	printf("Centroid1: %f %f %f \n",CC.centroidOfData1.x,CC.centroidOfData1.y,CC.centroidOfData1.z);
	printf("Centroid2: %f %f %f \n",CC.centroidOfData2.x,CC.centroidOfData2.y,CC.centroidOfData2.z);

	return CC;
}



gsl_matrix* CentroidCorrection(gsl_matrix *data, Dim3Point center, int size)
{
	int i;
	float holder;
	for (i = 0; i < size; ++i)
	{
		holder =  gsl_matrix_get (data, i, 0);
		holder -= center.x;
		gsl_matrix_set (data, i, 0, holder);

		holder =  gsl_matrix_get (data, i, 1);
		holder -= center.y;
		gsl_matrix_set (data, i, 1, holder);

		holder =  gsl_matrix_get (data, i, 2);
		holder -= center.z;
		gsl_matrix_set (data, i, 2, holder);

	}
	
	return data;
	

}


Transformation findRotationAndTranslation(gsl_matrix  *data1, gsl_matrix *data2, int size)
{
	Centroid CC;
	Transformation final;
	gsl_matrix *zeroMeanData1;
	gsl_matrix *zeroMeanData2;
	gsl_matrix *H,*V,*R,*c1,*c2;
	gsl_vector *work,*S;

	zeroMeanData1 = gsl_matrix_alloc (size, DIM);
	zeroMeanData2 = gsl_matrix_alloc (size, DIM);
	H = gsl_matrix_alloc(DIM,DIM);
	V = gsl_matrix_alloc(DIM,DIM);
	R = gsl_matrix_alloc(DIM,DIM);
	c1 = gsl_matrix_alloc(DIM,1);
	c2 = gsl_matrix_alloc(DIM,1);
	S = gsl_vector_alloc (DIM);
	work = gsl_vector_alloc (DIM);

	CC = findCentroid(data1, data2, size); // Find the centroid of the datasets.

	/* 
		Subtract the centroid from data points - i.e make the data set zero mean
	*/
	gsl_matrix_set (c1,0,0,CC.centroidOfData1.x);
	gsl_matrix_set (c1,1,0,CC.centroidOfData1.y);
	gsl_matrix_set (c1,2,0,CC.centroidOfData1.z);

	gsl_matrix_set (c2,0,0,CC.centroidOfData2.x);
	gsl_matrix_set (c2,1,0,CC.centroidOfData2.y);
	gsl_matrix_set (c2,2,0,CC.centroidOfData2.z);
	
	zeroMeanData1 = CentroidCorrection(data1, CC.centroidOfData1, size);
	zeroMeanData2 = CentroidCorrection(data2, CC.centroidOfData2, size);

	
	display(zeroMeanData1, size, DIM, "Zero Mean Data Point 1");
	display(zeroMeanData2, size, DIM, "Zero Mean Data Point 2");

	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, zeroMeanData1, zeroMeanData2, 0.0, H);


	display(H,DIM,DIM, "H");

	gsl_linalg_SV_decomp (H, V, S, work);

	display(H,DIM,DIM, "U");
	display(V,DIM,DIM, "V");

	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, V, H, 0.0, R);	

	//display(R,DIM,DIM, "Rotation Matrix");

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, R, c1, 1.0, c2);

	//display(c2,DIM,1,"Translation vector");

	final.rotationMatrix = R;
	final.translationVector = c2;

	return final;
}


