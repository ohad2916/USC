#include <string.h>
#include "spkmeans.h"
#pragma warning(disable:4996)


int main(int argc, char* argv[]) {
	char char_surfer = 'f';
	size_t commas = 0;
	size_t dimension = 1;
	size_t point_count = 0;
	size_t j,i;
	int return_value = 0;
	int mode;
	double** points_list;
	double* p;
	EIGEN_VALUES_VECTORS* jacobis_res;
	FILE* file;
	MATRIX weight_adj;MATRIX diag_degree_mat;MATRIX graph_laplacian;
	/*handeling cmdline arguments*/
	if (argc != 3) {
		printf("An Error Has Occured");
		return 1;
	}
	file = fopen(argv[2], "r");
	if (!file) {
		printf("an Error has occured!");
		return 1;
	}
	if (strcmp(argv[1], "wam") == 0)
		mode = 0;
	else if (strcmp(argv[1], "ddg") == 0)
		mode = 1;
	else if (strcmp(argv[1], "gl") == 0)
		mode = 2;
	else if (strcmp(argv[1], "jacobi") == 0)
		mode = 3;
	else {
		printf("an Error has occured!");
		return 1;
	}
	/*reading the file*/
	while ((char_surfer = fgetc(file)) != '\n' && char_surfer != EOF) {
		dimension += (char_surfer == ',');
	}
	commas = dimension - 1;
	while ((char_surfer = fgetc(file)) != EOF) {
		commas += (char_surfer == ',');
	}
	point_count = commas / (dimension - 1);
	/*printf("dimension is: %lu\nno.points:%lu\n", dimension, point_count);*/
	rewind(file);

	if (!(points_list = malloc(sizeof(double*) * point_count))) {
		printf("an Error has occured!");
		return 1;
	}
	if (!(p = malloc(sizeof(double) * point_count * dimension))) {
		printf("an Error has occured!");
		return 1;
	}
	for (i = 0; i < point_count; i++) {
		points_list[i] = p + (dimension * i);
	}

	for (i = 0; i < point_count; i++) {
		for (j = 0; j < dimension - 1; j++) {
			if (!fscanf(file, "%lf,", &points_list[i][j])) {
				freeMatrix(points_list);
				return 1;
			}
		}
		if (!fscanf(file, "%lf\n", &points_list[i][dimension - 1])) {
			freeMatrix(points_list);
			return 1;
		}
	}
	fclose(file);
	/*end reading file*/
	if (mode == 3) {
		jacobis_res = jacobis(points_list, point_count, 100, 0.00001);
		if (jacobis_res) {
			printVector(jacobis_res->eigen_values, point_count);
			printMatrix(jacobis_res->eigen_vectors_as_columns, point_count, point_count);
		}
		else {
			return_value = 1;
			printf("an Error has occured!");
		}
		freeMatrix(points_list);
		free(jacobis_res->eigen_values);
		freeMatrix(jacobis_res->eigen_vectors);
		freeMatrix(jacobis_res->eigen_vectors_as_columns);
		free(jacobis_res);
		return return_value;
	}
	weight_adj = weightedAdjacencyMatrix(points_list, point_count, dimension);
	if (mode == 0) {
		if (weight_adj)
			printMatrix(weight_adj, point_count, point_count);
		else {
			return_value = 1;
			printf("an Error has occured!");
		}
		freeMatrix(weight_adj);
		freeMatrix(points_list);
		return return_value;
	}
	diag_degree_mat = diagonalDegreeMatrix(weight_adj, point_count);
	if (mode == 1) {
		if (diag_degree_mat)
			printMatrix(diag_degree_mat, point_count, point_count);
		else {
			return_value = 1;
			printf("an Error has occured!");
		}
		freeMatrix(diag_degree_mat);
		freeMatrix(weight_adj);
		freeMatrix(points_list);
		return return_value;
	}
	graph_laplacian = graphLaplacian(weight_adj, diag_degree_mat, point_count);
	if (mode == 2) {
		if (graph_laplacian)
			printMatrix(graph_laplacian, point_count, point_count);
		else {
			return_value = 1;
			printf("an Error has occured!");
		}
		freeMatrix(graph_laplacian);
		freeMatrix(diag_degree_mat);
		freeMatrix(weight_adj);
		freeMatrix(points_list);
		return return_value;
	}

	/*cleanup*/
	freeMatrix(graph_laplacian);
	freeMatrix(diag_degree_mat);
	freeMatrix(weight_adj);
	freeMatrix(points_list);
	
	return return_value;



}