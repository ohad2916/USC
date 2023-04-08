#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MATRIX double**
#define VECTOR double*
#define _is_negative_zero(x) ((x == 0.0 && x < 0.0) || \
                              (x > -0.0001 && x < 0.0))

double twoNorm(double* a, double* b, size_t dimension);
typedef struct eigen_v {
	MATRIX eigen_vectors;
	MATRIX eigen_vectors_as_columns;
	VECTOR eigen_values;
}EIGEN_VALUES_VECTORS;
typedef struct k_eigenvs {
	MATRIX eigen_vectors;
	size_t k;
}FIRST_K_EIGEN;
typedef struct pivot_values_ {
	int i, j;
	double s, c;
}PIVOT_VALUES;

int freeMatrix(MATRIX a);
MATRIX allocateSquareMatrix(size_t n);
PIVOT_VALUES* generatePivotValues(MATRIX _A, size_t n);
int generatePivotMatrix(MATRIX target,size_t n,PIVOT_VALUES* pivot_values);
double vectorDot(VECTOR x, VECTOR y, size_t n);
MATRIX matrixDot(MATRIX a, MATRIX b, size_t n);
MATRIX transposeMatrix(MATRIX p,size_t n,const char* how);
EIGEN_VALUES_VECTORS* jacobis(MATRIX _A,size_t n,size_t iteration_limit, double epsilon);
size_t find_k(VECTOR e_values,size_t n);
int sortEigenVectors(EIGEN_VALUES_VECTORS* subject,size_t n);
MATRIX weightedAdjacencyMatrix(double** points, size_t n,size_t dimension);
MATRIX diagonalDegreeMatrix(double** adj_matrix, size_t n);
MATRIX graphLaplacian(double** adj_matrix, double** diag_matrix, size_t n);
MATRIX allocateNonSquareMatrix(size_t rows, size_t cols);
int printMatrix(double** mat, size_t rows, size_t cols) ;
int printVector(VECTOR v, size_t n);
