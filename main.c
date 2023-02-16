#include <stdio.h>
#include <stdlib.h>
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
	print_matrix(points_list, point_count, dimension);
	/*allocating weighted adjacency matrix*/
	double** adj_matrix;
	double* adj_data;
	if (!(adj_matrix = malloc(sizeof(double*) * point_count))) {
		printf("eror allocating");
		return 1;
	}
	if (!(adj_data = malloc(sizeof(double) * point_count * point_count))) {
		printf("eror allocating");
		return 1;
	}
	for (i = 0; i < point_count; i++) {
		adj_matrix[i] = adj_data + i * point_count;
	}
	/*calculating the matrix*/
	weighted_adj_calc(points_list, adj_matrix, point_count, dimension);
	print_matrix(adj_matrix, point_count, point_count);

	/*allocating diagonal degree matrix*/
	double** diag_matrix;
	double* diag_data;
	if (!(diag_matrix = malloc(sizeof(double*) * point_count))) {
		printf("eror allocating");
		return 1;
	}
	if (!(diag_data = calloc(point_count*point_count,sizeof(double)))){
		printf("eror allocating");
		return 1;
	}
	for (i = 0; i < point_count; i++) {
		diag_matrix[i] = diag_data + i * point_count;
	}
	/*calculating diagonal matrix*/
	diagonal_deg_calc(adj_matrix, diag_matrix, point_count);
	print_matrix(diag_matrix, point_count, point_count);

	/*allocating graph laplacian matrix*/
	double** laplac_matrix;
	double* laplac_data;
	if (!(laplac_matrix = malloc(sizeof(double*) * point_count))) {
		printf("eror allocating");
		return 1;
	}
	if (!(laplac_data = calloc(point_count * point_count, sizeof(double)))) {
		printf("eror allocating");
		return 1;
	}
	for (i = 0; i < point_count; i++) {
		laplac_matrix[i] = laplac_data + i * point_count;
	}
	/*calculating graph laplacian*/
	graph_laplacian_calc(adj_matrix, diag_matrix, laplac_matrix, point_count);
	print_matrix(laplac_matrix, point_count, point_count);


	return 1;
}