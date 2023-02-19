#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MATRIX double**
#define VECTOR double*

double twoNorm(double* a, double* b, size_t dimension) {
	size_t i = 0;
	double sum = 0;
	for (i = 0; i < dimension; i++) {
		sum += pow(a[i] - b[i], 2);
	}
	return sum;
}

typedef struct m_list{
	MATRIX this;
	struct m_list* next;
}MATRIX_LIST;
typedef struct eigen_v {
	MATRIX eigen_vectors;
	MATRIX eigen_vectors_as_columns;
	VECTOR eigen_values;
}EIGEN_VALUES_VECTORS;
typedef struct k_eigenvs {
	MATRIX eigen_vectors;
	size_t k;
}FIRST_K_EIGEN;

int append(MATRIX a, MATRIX_LIST* list) {
	if (!list->this)
		list->this = a;
	else {
		while (list->next)
			list = list->next;
		MATRIX_LIST* new = malloc(sizeof(MATRIX_LIST));
		new->next = NULL;
		new->this = a;
		list->next = new;
	}

}

MATRIX allocateSquareMatrix(size_t n) {
	double** matrix;
	double* data;
	size_t i;
	if (!(matrix = malloc(sizeof(double*) * n))) {
		printf("error allocating");
		return 1;
	}
	if (!(data = calloc(n*n, sizeof(double)))) {
		printf("error allocating");
		return 1;
	}
	for (i = 0; i <n; i++) {
		matrix[i] = data + i * n;
	}
	return matrix;
}

MATRIX generatePivotMatrix(MATRIX _A,size_t n) {
	MATRIX _P = allocateSquareMatrix(n);
	/*finding arg(s)max A(i,j)*/
	size_t i, j;
	size_t argmax_i = 0, argmax_j = 0;
	double max_value = -1.0;
	double phi,c,t,s;
	int sign = 0;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			if (fabs(_A[i][j]) > max_value) {
				argmax_i = i;
				argmax_j = j;
				max_value = fabs(_A[i][j]);
			}
		}
	}
	/*obtaining c,t*/
	phi = (_A[argmax_j][argmax_j] - _A[argmax_i][argmax_i]) / (2.0 * _A[argmax_i][argmax_j]);
	sign += (phi >= 0);
	sign -= (sign == 0);
	t = sign / ((abs(phi) + sqrt(phi * phi + 1)));
	c = 1 / (sqrt(t * t + 1));
	s = t * c;
	/*ammending pivot matrix*/
	for (i = 0; i < n; i++) {
		_P[i][i] = 1;
	}
	_P[argmax_i][argmax_i] = c;
	_P[argmax_i][argmax_j] = s;
	_P[argmax_j][argmax_j] = c;
	if(s != 0.0)
		_P[argmax_j][argmax_i] = -s;
	
	return _P;
}

double vectorDot(VECTOR x, VECTOR y, size_t n) {
	size_t i;
	double sum = 0;
	for (i = 0; i < n; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}

MATRIX matrixDot(MATRIX a, MATRIX b, size_t n) {
	size_t i, j,k;
	MATRIX dot = allocateSquareMatrix(n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0;k < n;k++) {
				dot[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return dot;
}

MATRIX transposeMatrix(MATRIX p,size_t n,const char* how){
	size_t i, j;
	double temp1;
	double temp2;
	MATRIX res;
	if (how != "in-place") {
		res = allocateSquareMatrix(n);
		for (i = 0;i < n;i++)
			res[i][i] = p[i][i];
	}
	else {
		res = p;
	}
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			temp1 = p[i][j];
			temp2 = p[j][i];
			res[i][j] = temp2;
			res[j][i] = temp1;
		}
	}
	return res;
}

EIGEN_VALUES_VECTORS* jacobis(MATRIX _A,size_t n,size_t iteration_limit, double epsilon) {
	EIGEN_VALUES_VECTORS* res = malloc(sizeof(EIGEN_VALUES_VECTORS));
	if (!res) {
		printf("an Error has Occured");
		return 1;
	}
	double current_off_diag_sum = -1.0;
	double new_off_diag_sum = .0;
	double convergence_diff = .0;
	size_t i, j,k;
	MATRIX dotted_pivots = allocateSquareMatrix(n);
	for (i = 0; i < n; i++) {
		dotted_pivots[i][i] = 1;
	}

	for (i = 0; i < iteration_limit; i++) {
		MATRIX pivot_matrix = generatePivotMatrix(_A, n);
		/*generating eigen vector matrix */
		MATRIX temp_pivot_pointer = matrixDot(dotted_pivots, pivot_matrix,n);
		free(*dotted_pivots);
		free(dotted_pivots);
		dotted_pivots = temp_pivot_pointer;
		MATRIX dot_res = matrixDot(_A, pivot_matrix, n);
		MATRIX transposed_pivot = transposeMatrix(pivot_matrix, n,"in-place");
		MATRIX second_dot_res = matrixDot(transposed_pivot, dot_res, n);
		if (i != 0) {
			free(*_A);
			free(_A);
		}
		_A = second_dot_res;
		//checking convergence
		new_off_diag_sum = 0.0;
		for (k = 0; k < n; k++) {
			for (j = n+1; j < n; j++) {
				new_off_diag_sum += _A[i][j]* _A[i][j];
			}
		}
		new_off_diag_sum *= 2;
		/*if (current_off_diag_sum - new_off_diag_sum <= epsilon)
			break;*/
		//cleanup
		free(*transposed_pivot);
		free(transposed_pivot);
		free(*dot_res);
		free(dot_res);
	}
	
	res->eigen_vectors_as_columns = dotted_pivots;
	res->eigen_vectors = transposeMatrix(dotted_pivots, n, "new");
	VECTOR eigen_values = malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		eigen_values[i] = _A[i][i];
	}
	res->eigen_values = eigen_values;
	free(*_A);
	free(_A);
	return res;
}

size_t find_k(VECTOR e_values,size_t n) {
	size_t i;
	double max_gap = .0;
	size_t argMax = 0;
	for ( i = 0; i < n - 1; i++) {
		if (fabs(e_values[i + 1] - e_values[i]) > max_gap) {
			max_gap = fabs(e_values[i + 1] - e_values[i]);
			argMax = i;
		}
	}
	return i;
}

int sortEigenVectors(EIGEN_VALUES_VECTORS* subject,size_t n) {
	size_t i, j;
	double temp_value;
	VECTOR temp_vector;
	MATRIX e_vectors = subject->eigen_vectors;
	VECTOR e_values = subject->eigen_values;
	for (i = 0; i < n - 1;i++) {
		for (j = 0; j < n - i -1; j++) {
			if (e_values[j] > e_values[j + 1]) {
				temp_value = e_values[j];
				temp_vector = e_vectors[j];
				e_values[j] = e_values[j + 1];
				e_vectors[j] = e_vectors[j + 1];
				e_values[j + 1] = temp_value;
				e_vectors[j + 1] = temp_vector;
			}
		}
	}
	free(*(subject->eigen_vectors_as_columns));
	free(subject->eigen_vectors_as_columns);
	subject->eigen_vectors_as_columns = transposeMatrix(e_vectors, n, "new");
	return 0;
}

MATRIX weightedAdjacencyMatrix(double** points, size_t n,size_t dimension) {
	MATRIX res = allocateSquareMatrix(n);
	size_t i;
	size_t j;
	for (i = 0; i < n; i++) {
		for (j = i; j < n;j++) {
			res[i][j] = twoNorm(points[i], points[j], dimension);
			res[j][i] = res[i][j];
		}
	}
	return res;
}

MATRIX diagonalDegreeMatrix(double** adj_matrix, size_t n) {
	MATRIX res = allocateSquareMatrix(n);
	size_t i;
	size_t j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			res[i][i] += adj_matrix[i][j];
		}
	}
	return res;
}

MATRIX graphLaplacian(double** adj_matrix, double** diag_matrix, size_t n) {
	MATRIX res = allocateSquareMatrix(n);
	size_t i, j;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			res[i][j] = diag_matrix[i][j] - adj_matrix[i][j];
			res[j][i] = res[i][j];
		}
	}
	return res;
}

int printMatrix(double** mat, size_t rows, size_t cols) {
	int i = 0; int j = 0;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			printf("%.2f", mat[i][j]);
			if (j < cols - 1)
				printf(",");
		}
		printf("\n");
	}
	return 0;
}
int printVector(VECTOR v, size_t n) {
	size_t i = 0;
	for (i = 0; i < n - 1;i++) {
		printf("%.2f,", v[i]);
	}
	printf("%.2f",v[n - 1]);
	printf("\n");
	return 0;
}