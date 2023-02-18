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
typedef struct eigen_v {
	MATRIX eigen_vectors;
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

EIGEN_VALUES_VECTORS* jacobis(MATRIX _A,size_t n,size_t iteration_limit, double epsilon) {
	EIGEN_VALUES_VECTORS* res = malloc(sizeof(EIGEN_VALUES_VECTORS));
	MATRIX_LIST p_list;
	p_list.this = NULL;
	p_list.next = NULL;
	double current_off_diag_sum = -1.0;
	double new_off_diag_sum = .0;
	double convergence_diff = .0;
	size_t i, j;

	for (i = 0; i < iteration_limit; i++) {
		MATRIX pivot_matrix = generate_pivot_matrix(_A, n);
		append(pivot_matrix, &p_list);
		MATRIX dot_res = matrix_dot(_A, pivot_matrix, n);
		MATRIX transposed_pivot = transpose_matrix(pivot_matrix, n);
		MATRIX second_dot_res = matrix_dot(transposed_pivot, dot_res, n);
		if(i != 0)
			free(*_A);
		_A = second_dot_res;
		new_off_diag_sum = 0.0;
		for (i = 0; i < n; i++) {
			for (j = n+1; j < n; j++) {
				new_off_diag_sum += _A[i][j]* _A[i][j];
			}
		}
		new_off_diag_sum *= 2;
		if (current_off_diag_sum - new_off_diag_sum <= epsilon)
			i = 100;
		free(*transposed_pivot);
		free(*dot_res);
	}
	MATRIX dot_res = allocate_square_matrix(n);
	for (i = 0; i < n; i++) {
		dot_res[i][i] = 1;
	}
	MATRIX temp_res = NULL;
	do {
		temp_res = matrix_dot(dot_res, p_list.this, 3);
		free(*dot_res);
		dot_res = temp_res;
		free(*p_list.this);
		if(p_list.next != NULL)
			p_list = *p_list.next;

	} while (p_list.next != NULL);
	res->eigen_vectors = transpose_matrix(dot_res,n);
	free(*dot_res);
	VECTOR eigen_values = malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		eigen_values[i] = _A[i][i];
	}
	res->eigen_values = eigen_values;
	return res;
}

int find_k(VECTOR e_values,size_t n) {
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

int sort_eigen_vectors(VECTOR* e_vectors, VECTOR e_values,size_t n) {
	size_t i, j;
	double temp_value;
	VECTOR temp_vector;
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
	return 0;
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
int print_vector(VECTOR v, size_t n) {
	size_t i = 0;
	for (i = 0; i < n;i++) {
		printf("%.2f,", v[i]);
	}
	printf("\n");
	return 0;
}