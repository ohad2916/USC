#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MATRIX double**
#define VECTOR double*
#define _is_negative_zero(x) ((x == 0.0 && x < 0.0) || \
                              (x > -0.0001 && x < 0.0))

double twoNorm(double* a, double* b, size_t dimension) {
	size_t i = 0;
	double sum = 0;
	for (i = 0; i < dimension; i++) {
		sum += pow(a[i] - b[i], 2);
	}
	return sum;
}

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
int freeMatrix(MATRIX a) {
	free(*a);
	free(a);
	return 0;
}


MATRIX allocateSquareMatrix(size_t n) {
	double** matrix;
	double* data;
	size_t i;
	if (!(matrix = malloc(sizeof(double*) * n))) {
		return NULL;
	}
	if (!(data = calloc(n*n, sizeof(double)))) {
		free(matrix);
		return NULL;
	}
	for (i = 0; i <n; i++) {
		matrix[i] = data + (i * n);
	}
	return matrix;
}

PIVOT_VALUES* generatePivotValues(MATRIX _A, size_t n) {	
	size_t i, j, argmax_i, argmax_j;
	double max_value = -1.0;
	double phi, c, t, s;
	int sign;
	PIVOT_VALUES* res = malloc(sizeof(PIVOT_VALUES));
	if (!res) {
		return NULL;
	}
	/*finding arg(s)max A(i,j)*/
	argmax_i = 0; argmax_j = 0;
	sign = 0;
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
	t = sign / ((fabs(phi) + sqrt(phi * phi + 1)));
	c = 1 / (sqrt(t * t + 1));
	s = t * c;
	res->i = argmax_i;res->j = argmax_j;res->s = s;res->c = c;
	return res;
}

int generatePivotMatrix(MATRIX target,size_t n,PIVOT_VALUES* pivot_values) {
	size_t i, j;
	double s, c;
	for (i = 0; i < n; i++) {
		target[i][i] = 1;
	}
	i = pivot_values->i;
	j = pivot_values->j;
	s = pivot_values->s;
	c = pivot_values->c;
	target[i][i] = c;
	target[i][j] = s;
	target[j][j] = c;
	if(s != 0.0)
		target[j][i] = -s;
	return 0;
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
	if (!dot)
		return NULL;
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
	if (strcmp(how,"in-place") != 0) {
		res = allocateSquareMatrix(n);
		if (!res)
			return NULL;
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
	double current_off_diag_sum = 2147483647;
	double new_off_diag_sum = .0;
	MATRIX dotted_pivots = NULL; MATRIX pivot_matrix = NULL;
	VECTOR eigen_values = NULL;VECTOR I_col;VECTOR J_col;
	double a_ii_tag, a_jj_tag,s,c;
	size_t i, j,k, _i,_j,m;
	PIVOT_VALUES* pivot_values = NULL;
	EIGEN_VALUES_VECTORS* res = malloc(sizeof(EIGEN_VALUES_VECTORS));
	if (!res) {
		return NULL;
	}
	pivot_matrix = allocateSquareMatrix(n);
	dotted_pivots = allocateSquareMatrix(n);
	if (!dotted_pivots) {
		free(res);
		return NULL;
	}
	for (i = 0; i < n; i++) {
		dotted_pivots[i][i] = 1;
	}
	if (!(J_col = malloc(sizeof(double) * n))) {
		free(res);
		freeMatrix(dotted_pivots);
		return NULL;
	}
	if (!(I_col = malloc(sizeof(double) * n))) {
		freeMatrix(dotted_pivots);
		free(J_col);
		free(res);
		return NULL;
	}

	for (i = 0; i < iteration_limit; i++) {
		if (!(pivot_values = generatePivotValues(_A, n))) {
			free(J_col);
			free(I_col);
			freeMatrix(dotted_pivots);
			free(res);
			return NULL;
		}

		generatePivotMatrix(pivot_matrix, n,pivot_values);

		_i = pivot_values->i;
		_j = pivot_values->j;
		s = pivot_values->s;
		c = pivot_values->c;
		/*generating eigen vector matrix */
		memset(I_col, 0.0, n*(sizeof(double)));
		memset(J_col, 0.0, n*(sizeof(double)));

		for (k = 0; k < n;k++) {
			for (m = 0; m < n; m++) {
				I_col[k] += dotted_pivots[k][m] * pivot_matrix[m][_i];
				J_col[k] += dotted_pivots[k][m] * pivot_matrix[m][_j];
			}
		}
		for (k = 0; k < n;k++) {
			dotted_pivots[k][_i] = I_col[k];
			dotted_pivots[k][_j] = J_col[k];
		}

		/*------calculating A'------*/
		
			/*copying columns*/
	
		for (k = 0; k < n;k++) {
			I_col[k] = _A[k][_i];
			J_col[k] = _A[k][_j];
		}
		a_ii_tag = pow(c, 2) * _A[_i][_i] + pow(s, 2) * _A[_j][_j] - 2 * s * c * _A[_i][_j];
		a_jj_tag = pow(s, 2) * _A[_i][_i] + pow(c, 2) * _A[_j][_j] + 2 * s * c * _A[_i][_j];
		for (k = 0; k < n;k++) {
			_A[k][ _i] = c * I_col[k] - s * J_col[k];
			_A[_i][k] = c * I_col[k] - s * J_col[k];
			_A[k][_j] = c * J_col[k] + s * I_col[k];
			_A[_j][k] = c * J_col[k] + s * I_col[k];
		}
		_A[_i][_i] = a_ii_tag;
		_A[_j][_j] = a_jj_tag;
		_A[_i][_j] = 0.0;
		_A[_j][_i] = 0.0;
		/*checking convergence*/
		new_off_diag_sum = 0.0;
		for (k = 0; k < n; k++) {
			for (j = k+1; j < n; j++) {
				new_off_diag_sum += _A[k][j] * _A[k][j];
			}
		}
		new_off_diag_sum *= 2.0;
		/*cleanup*/
		/*free(*transposed_pivot);
		free(transposed_pivot);
		free(*dot_res);
		free(dot_res);*/

		/*zeroing pivot matrix*/
		pivot_matrix[_i][_i] = 0.0;
		pivot_matrix[_i][_j] = 0.0;
		pivot_matrix[_j][_i] = 0.0;
		pivot_matrix[_j][_j] = 0.0;
		/*breaking in case converged*/
		if ((current_off_diag_sum - new_off_diag_sum) <= epsilon) {
			break;
		}
		current_off_diag_sum = new_off_diag_sum;
		free(pivot_values);
		pivot_values = NULL;
	}
	
	res->eigen_vectors_as_columns = dotted_pivots;
	eigen_values = malloc(sizeof(double) * n);
	if (!eigen_values) {
		free(J_col);
		free(I_col);
		freeMatrix(dotted_pivots);
		free(res);
		return NULL;
	}
	res->eigen_vectors = transposeMatrix(dotted_pivots, n, "new");
	if (!(res->eigen_vectors)) {
		freeMatrix(dotted_pivots);
		free(eigen_values);
		free(res);
		return NULL;
	}
	for (i = 0; i < n; i++) {
		eigen_values[i] = _A[i][i];
	}
	for (i = 0; i < n;i++) {
		if (_is_negative_zero(eigen_values[i])) {
			for (j = 0;j < n;j++) {
				dotted_pivots[j][i] *= -1.0;
			}
			eigen_values[i] = 0.0;
		}
	}
	res->eigen_values = eigen_values;

	free(J_col);
	free(I_col);
	free(pivot_values);
	freeMatrix(pivot_matrix);
	return res;
}

size_t find_k(VECTOR e_values,size_t n) {
	size_t i;
	double max_gap = .0;
	size_t argMax = 0;
	for ( i = 0; i < n/2; i++) {
		if (e_values[i + 1] - e_values[i] > max_gap) {
			max_gap = e_values[i + 1] - e_values[i];
			argMax = i;
		}
	}
	return argMax + 1;
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
	size_t i;
	size_t j;
	MATRIX res = allocateSquareMatrix(n);
	if (!res) {
		return NULL;
	}
	for (i = 0; i < n;i++) {
		res[i][i] = 0;
	}
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n;j++) {
			double val = twoNorm(points[i], points[j], dimension);
			val *= -1;
			val /= 2;
			val = exp(val);
			res[i][j] = val;
			res[j][i] = val;
		}
	}
	return res;
}

MATRIX diagonalDegreeMatrix(double** adj_matrix, size_t n) {
	size_t i;
	size_t j;
	MATRIX res = allocateSquareMatrix(n);
	if (!res)
		return NULL;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			res[i][i] += adj_matrix[i][j];
		}
	}
	return res;
}

MATRIX graphLaplacian(double** adj_matrix, double** diag_matrix, size_t n) {
	size_t i, j;
	MATRIX res = allocateSquareMatrix(n);
	if (!res)
		return NULL;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			res[i][j] = diag_matrix[i][j] - adj_matrix[i][j];
			res[j][i] = res[i][j];
		}
	}
	return res;
}

MATRIX allocateNonSquareMatrix(size_t rows, size_t cols) {
	MATRIX matrix;
	double* data;
	size_t i;
	if (!(matrix = malloc(sizeof(double*) * rows))) {
		return NULL;
	}
	if (!(data = calloc(rows*cols, sizeof(double)))) {
		free(matrix);
		return NULL;
	}
	for (i = 0; i < rows; i++) {
		matrix[i] = data + i * cols;
	}
	return matrix;
}

int printMatrix(double** mat, size_t rows, size_t cols) {
	size_t i = 0; size_t j = 0;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols - 1; j++) {
			printf("%.4f,", mat[i][j]);
		}
		printf("%.4f\n", mat[i][cols - 1]);
	}
	return 0;
}
int printVector(VECTOR v, size_t n) {
	size_t i = 0;
	for (i = 0; i < n - 1;i++) {
		printf("%.4f,", v[i]);
	}
	printf("%.4f",v[n - 1]);
	printf("\n");
	return 0;
}

