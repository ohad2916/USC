#include <math.h>
#include <stdio.h>

double two_norm(double* a, double* b, size_t dimension) {
	size_t i = 0;
	double sum = 0;
	for (i = 0; i < dimension; i++) {
		sum += pow(a[i] - b[i], 2);
	}
	return sum;
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