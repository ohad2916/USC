#include <math.h>
#include <stdio.h>
#define MATRIX double**
#define VECTOR double*

double two_norm(double* a, double* b, size_t dimension) {
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

MATRIX allocate_square_matrix(size_t n) {
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

MATRIX generate_pivot_matrix(MATRIX _A,size_t n) {
	MATRIX _P = allocate_square_matrix(n);
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

double vector_dot(VECTOR x, VECTOR y, size_t n) {
	size_t i;
	double sum = 0;
	for (i = 0; i < n; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}

MATRIX matrix_dot(MATRIX a, MATRIX b, size_t n) {
	size_t i, j,k;
	MATRIX dot = allocate_square_matrix(n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0;k < n;k++) {
				dot[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return dot;
}

MATRIX transpose_matrix(MATRIX p,size_t n) {
	size_t i, j;
	double temp;
	MATRIX res = allocate_square_matrix(n);
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			res[j][i] = p[i][j];
			res[i][j] = p[j][i];
		}
		res[i][i] = p[i][i];
	}
	return res;
}


int weighted_adj_calc(double** points, double** target, size_t no_points,size_t dimension) {
	size_t i;
	size_t j;
	for (i = 0; i < no_points; i++) {
		for (j = i; j < no_points;j++) {
			target[i][j] = two_norm(points[i], points[j], dimension);
			target[j][i] = target[i][j];
		}
	}
	return 0;
}

int diagonal_deg_calc(double** adj_matrix, double** target, size_t no_points) {
	size_t i;
	size_t j;
	for (i = 0; i < no_points; i++) {
		for (j = 0; j < no_points; j++) {
			target[i][i] += adj_matrix[i][j];
		}
	}
	return 0;
}

int graph_laplacian_calc(double** adj_matrix, double** diag_matrix, double** target, size_t no_points) {
	size_t i, j;
	for (i = 0; i < no_points; i++) {
		for (j = i; j < no_points; j++) {
			target[i][j] = diag_matrix[i][j] - adj_matrix[i][j];
			target[j][i] = target[i][j];
		}
	}
}

int print_matrix(double** mat, size_t rows, size_t cols) {
	int i = 0; int j = 0;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			printf("%.2f,", mat[i][j]);
		}
		printf("\n");
	}
	return 0;
}