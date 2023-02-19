#include <string.h>
#include "spkmeans.h"
#pragma warning(disable:4996)

int main(int argc, char* argv[]) {
	//handeling cmdline arguments
	FILE* file = fopen(argv[2], "r");
	if (!file) {
		printf("an Error has occured!");
		return 1;
	}
	int mode;
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
	char char_surfer = 'f';
	size_t commas = 0;
	size_t dimension = 1;
	size_t point_count = 0;
	size_t j;
	size_t i;
	while ((char_surfer = fgetc(file)) != '\n' && char_surfer != EOF) {
		dimension += (char_surfer == ',');
	}
	commas = dimension - 1;
	while ((char_surfer = fgetc(file)) != EOF) {
		commas += (char_surfer == ',');
	}
	point_count = commas / (dimension - 1);
	printf("dimension is: %d\nno.points:%d\n", dimension, point_count);
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
		points_list[i] = p + (dimension * i);
	}

	for (i = 0; i < point_count; i++) {
		for (j = 0; j < dimension - 1; j++) {
			if (!fscanf(file, "%lf,", &points_list[i][j]))
				return 1;
		}
		if (!fscanf(file, "%lf\n", &points_list[i][dimension - 1]))
			return 1;
	}
	/*end reading file*/
	MATRIX weight_adj = weightedAdjacencyMatrix(points_list, point_count, dimension);
	if (mode == 0) {
		printf("GOAL 0\n");
		printMatrix(weight_adj, point_count, point_count);
		free(*weight_adj);
		free(weight_adj);
		free(*points_list);
		free(points_list);
		return 0;
	}
	MATRIX diag_degree_mat = diagonalDegreeMatrix(weight_adj, point_count);
	if (mode == 1) {
		printf("GOAL 1\n");
		printMatrix(diag_degree_mat, point_count, point_count);
		free(*diag_degree_mat);
		free(diag_degree_mat);
		free(*weight_adj);
		free(weight_adj);
		free(*points_list);
		free(points_list);
		return 0;
	}
	MATRIX graph_laplacian = graphLaplacian(weight_adj, diag_degree_mat, point_count);
	if (mode == 2) {
		printf("GOAL 2\n");
		printMatrix(graph_laplacian,point_count,point_count);
		free(*graph_laplacian);
		free(graph_laplacian);
		free(*diag_degree_mat);
		free(diag_degree_mat);
		free(*weight_adj);
		free(weight_adj);
		free(*points_list);
		free(points_list);
		return 0;
	}
	printf("GOAL 3\n");
	EIGEN_VALUES_VECTORS* jacobis_res = jacobis(graph_laplacian, point_count, 100, 0.000001);
	printVector(jacobis_res->eigen_values,point_count);
	printMatrix(jacobis_res->eigen_vectors_as_columns, point_count, point_count);

	/*cleanup*/
	free(*graph_laplacian);
	free(graph_laplacian);
	free(*diag_degree_mat);
	free(diag_degree_mat);
	free(*weight_adj);
	free(weight_adj);
	free(*points_list);
	free(points_list);
	free(jacobis_res->eigen_values);
	free(*(jacobis_res->eigen_vectors));
	free(jacobis_res->eigen_vectors);
	free(*(jacobis_res->eigen_vectors_as_columns));
	free(jacobis_res->eigen_vectors_as_columns);
	free(jacobis_res);
	return 0;



}