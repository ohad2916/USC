#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spkmeans.h"
#pragma warning(disable:4996)

int main() {
	FILE* file = fopen("Test Files/input_0.txt","r");
	if (!file) {
		printf("error opening file");
		return 1;
	}
	char char_surfer = 'f';
	size_t dimension = 1;
	size_t point_count = 1;
	size_t j;
	size_t i;
	while ((char_surfer = fgetc(file)) != '\n') {
		dimension += (char_surfer == ',');
	}
	while ((char_surfer = fgetc(file)) != EOF) {
		point_count += (char_surfer == '\n');
	}
	printf("dimension is: %d\nno.points:%d\n", dimension,point_count);
	rewind(file);

	double** points_list;
	if (!(points_list = malloc(sizeof(double*) * point_count))) {
		printf("error allocating");
		return 1;
	}
	double* p;
	if (!(p = malloc(sizeof(double) * point_count * dimension))) {
		printf("error allocating");
		return 1;
	}
	for (i = 0; i < point_count; i++) {
		points_list[i] = p + dimension * i;
	}
	
	for (i = 0; i < point_count; i++) {
		for (j = 0; j < dimension -1; j++) {
			if (!fscanf(file, "%lf,", &points_list[i][j]))
				return 1;
		}
		if (!fscanf(file, "%lf\n", &points_list[i][dimension - 1]))
			return 1;
	}
	//print_matrix(points_list, point_count, dimension);
	/*allocating weighted adjacency matrix*/
	MATRIX adj_matrix = allocate_square_matrix(point_count);
	/*calculating the matrix*/
	weighted_adj_calc(points_list, adj_matrix, point_count, dimension);
	//print_matrix(adj_matrix, point_count, point_count);

	/*allocating diagonal degree matrix*/
	MATRIX diag_matrix = allocate_square_matrix(point_count);
	/*calculating diagonal matrix*/
	diagonal_deg_calc(adj_matrix, diag_matrix, point_count);
	//print_matrix(diag_matrix, point_count, point_count);

	/*allocating graph laplacian matrix*/
	MATRIX laplac_matrix = allocate_square_matrix(point_count);
	/*calculating graph laplacian*/
	graph_laplacian_calc(adj_matrix, diag_matrix, laplac_matrix, point_count);
	print_matrix(laplac_matrix, point_count, point_count);

	/*debuging pivot*/
	MATRIX pivot_matrix = generate_pivot_matrix(laplac_matrix, point_count);
	print_matrix(pivot_matrix, point_count, point_count);
	printf("\n\n");
	/*jacobis*/
	MATRIX _A = allocate_square_matrix(3);
	_A[0][0] = 1;_A[0][1] = 0;_A[0][2] = 0;
	_A[1][0] = 0;_A[1][1] = 4;_A[1][2] = -1;
	_A[2][0] = 0;_A[2][1] = -1;_A[2][2] = 2;

	EIGEN_VALUES_VECTORS* test_jacobi = jacobis(_A, 3, 100, 0.0);

	print_matrix(test_jacobi->eigen_vectors, 3, 3);
	printf("\n");
	print_vector(test_jacobi->eigen_values, 3);
	printf("\n");
	sort_eigen_vectors(test_jacobi->eigen_vectors, test_jacobi->eigen_values, 3);
	print_matrix(test_jacobi->eigen_vectors, 3, 3);
	printf("\n");
	print_vector(test_jacobi->eigen_values,3);


	return 1;
}