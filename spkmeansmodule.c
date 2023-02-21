#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
#define PYTHON_MODE


double euc_d(double* p, double* q, size_t dim) {
    double d_sqrd_sum = 0;
    size_t i = 0;
    while (i < dim) {
        d_sqrd_sum += pow(p[i] - q[i], 2);
        i++;
    }
    return sqrt(d_sqrd_sum);
}
int free_memory(double* a, double* b, double* c, double** data, double** mean, double** curr) {
    free(a);
    free(b);
    free(c);
    free(data);
    free(mean);
    free(curr);
    return 0;
}

static PyObject* fit(PyObject* self, PyObject* args) {
    PyObject* data_lst;
    PyObject* centroid_lst;
    PyObject* point;
    size_t iter;
    double epsilon;
    if (!PyArg_ParseTuple(args, "OOnd", &data_lst, &centroid_lst, &iter, &epsilon)) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    /*data_lst = dataframe */
    size_t no_points = PyObject_Length(data_lst);
    point = PyList_GetItem(data_lst, 0);
    size_t dimension = PyObject_Length(point);
    size_t centroid_size = PyObject_Length(centroid_lst);

    double* p = NULL;
    double** data = NULL;
    double* b = NULL;
    double** cluster_mean = NULL;
    double* c = NULL;
    double** new_cluster = NULL;
    size_t j, m, i;
    data = calloc(no_points, sizeof(double*));
    if (!data) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        free_memory(b, c, p, new_cluster, cluster_mean, data);
        return NULL;
    }
    p = NULL;
    p = (double*)calloc(no_points * dimension, sizeof(double));
    if (!p) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        free_memory(b, c, p, new_cluster, cluster_mean, data);
        return NULL;
    }
    for (i = 0; i < no_points; i++) {
        data[i] = p + i * dimension;
    }

    for (i = 0; i < no_points; i++) {
        point = PyList_GetItem(data_lst, i);
        for (j = 0; j < dimension; j++) {
            data[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
            //printf("%f,", data[i][j]);
        }
    }
    /*convert current cluster list to a c array*/
    b = calloc(centroid_size * dimension, sizeof(double));
    cluster_mean = calloc(centroid_size, sizeof(double*));
    if (!b || !cluster_mean) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        free_memory(b, c, p, new_cluster, cluster_mean, data);
        return NULL;
    }
    for (i = 0; i < centroid_size; i++) {
        cluster_mean[i] = b + i * (dimension);
    }

    for (i = 0; i < centroid_size; i++) {
        point = PyList_GetItem(centroid_lst, i);
        for (j = 0; j < dimension; j++) {
            cluster_mean[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    /*main algorithm*/
    double** curr_X = data;
    /*allocate temporary clusters to decide convergence*/
    c = calloc(centroid_size * (dimension + 1), sizeof(double));
    new_cluster = calloc(centroid_size, sizeof(double*));
    if (!c || !new_cluster) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        free_memory(b, c, p, new_cluster, cluster_mean, data);
        return NULL;
    }
    for (i = 0; i < centroid_size; i++) {
        new_cluster[i] = c + i * (dimension + 1);
    }

    i = 0;
    int converged = 0;
    while (i < iter && !converged) {
        /*zero out new cluster array*/
        memset(c, 0, centroid_size * (dimension + 1) * sizeof(double));
        /*decide closest cluster against the original and update it with new Xi*/
        for (m = 0; m < no_points; m++) {
            size_t min_cluster_index = 0;
            double min_value = INT32_MAX;
            for (j = 0; j < centroid_size; j++) {
                double curr_euc_d = euc_d(cluster_mean[j], *(curr_X + m), dimension);
                if (curr_euc_d < min_value) {
                    min_value = curr_euc_d;
                    min_cluster_index = j;
                }
            }
            /*updating new cluster, just adding for now. divide later.*/
            double* min_cluster = new_cluster[min_cluster_index];
            for (j = 0; j < dimension; j++) {
                min_cluster[j] += curr_X[m][j];
            }
            min_cluster[dimension]++;
        }
        /*calculate the actual means*/
        for (j = 0; j < centroid_size; j++) {
            for (m = 0; m < dimension; m++) {
                new_cluster[j][m] /= new_cluster[j][dimension];
            }
        }

        /*decide convegerence*/
        double max_Duk = 0;
        double curr_Muk = 0;
        for (j = 0; j < centroid_size; j++) {
            curr_Muk = euc_d(cluster_mean[j], new_cluster[j], dimension);
            if (curr_Muk > max_Duk)
                max_Duk = curr_Muk;
        }
        if (max_Duk <= epsilon) {
            /*print statement for debugging*/
            /*printf("Converged after: %d iterations\n", (int)i + 1);*/
            converged = 1;
        }
        i++;
        /*copy new cluster to old ones*/
        for (j = 0; j < centroid_size; j++) {
            for (m = 0; m < dimension; m++) {
                cluster_mean[j][m] = new_cluster[j][m];
            }
        }
    }
    PyObject* py_centroids = PyList_New(centroid_size);
    for (m = 0; m < centroid_size; m++) {
        PyObject* m_cluster = PyList_New(dimension);
        PyList_SetItem(py_centroids, m, m_cluster);        //raises exception
        for (j = 0; j < dimension; j++) {
            PyList_SetItem(m_cluster, j, PyFloat_FromDouble(cluster_mean[m][j]));
        }
    }
    free_memory(b, c, p, new_cluster, cluster_mean, data);
    return py_centroids;
}
static PyObject* wam(PyObject* self, PyObject* args) {
    PyObject* py_point_list;
    if (!PyArg_ParseTuple(args, "O", &py_point_list)) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    size_t no_points = PyObject_Length(py_point_list);
    PyObject* point = PyList_GetItem(py_point_list, 0);
    size_t dimension = PyObject_Length(point);
    size_t i, j;
    //allocate c matrix
    MATRIX c_point_list = allocateNonSquareMatrix(no_points, dimension);
    if (!c_point_list) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < no_points; i++) {
        point = PyList_GetItem(py_point_list, i);
        for (j = 0; j < dimension; j++) {
            c_point_list[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    MATRIX weighted_adj_matrix = weightedAdjacencyMatrix(c_point_list, no_points, dimension);
    if (!weighted_adj_matrix) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);
        return NULL;
    }
    PyObject* py_wam = PyList_New(no_points);
    for (i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(py_wam, i, row);        //raises exception
        for (j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(weighted_adj_matrix[i][j]));
        }
    }
    freeMatrix(c_point_list);
    freeMatrix(weighted_adj_matrix);
    return py_wam;
}
static PyObject* ddg(PyObject* self, PyObject* args) {
    PyObject* py_point_list;
    if (!PyArg_ParseTuple(args, "O", &py_point_list)) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    size_t no_points = PyObject_Length(py_point_list);
    PyObject* point = PyList_GetItem(py_point_list, 0);
    size_t dimension = PyObject_Length(point);
    size_t i, j;
    //allocate c matrix
    MATRIX c_point_list = allocateNonSquareMatrix(no_points, dimension);
    if (!c_point_list) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < no_points; i++) {
        point = PyList_GetItem(py_point_list, i);
        for (j = 0; j < dimension; j++) {
            c_point_list[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    MATRIX weighted_adj_matrix = weightedAdjacencyMatrix(c_point_list, no_points, dimension);
    if (!weighted_adj_matrix) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);
        return NULL;
    }
    MATRIX diagonal_deg_matrix = diagonalDegreeMatrix(weighted_adj_matrix,no_points);
    if (!diagonal_deg_matrix) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);
        freeMatrix(weighted_adj_matrix);
        return NULL;
    }
    PyObject* py_deg = PyList_New(no_points);
    for (i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(py_deg, i, row);        //raises exception
        for (j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(diagonal_deg_matrix[i][j]));
        }
    }
    freeMatrix(c_point_list);
    freeMatrix(weighted_adj_matrix);
    freeMatrix(diagonal_deg_matrix);
    return py_deg;
}
static PyObject* gl(PyObject* self, PyObject* args) {
    PyObject* py_point_list;
    if (!PyArg_ParseTuple(args, "O", &py_point_list)) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    size_t no_points = PyObject_Length(py_point_list);
    PyObject* point = PyList_GetItem(py_point_list, 0);
    size_t dimension = PyObject_Length(point);
    size_t i, j;
    //allocate c matrix
    MATRIX c_point_list = allocateNonSquareMatrix(no_points, dimension);
    if (!c_point_list) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < no_points; i++) {
        point = PyList_GetItem(py_point_list, i);
        for (j = 0; j < dimension; j++) {
            c_point_list[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    MATRIX weighted_adj_matrix = weightedAdjacencyMatrix(c_point_list, no_points, dimension);
    if (!weighted_adj_matrix) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);

        return NULL;
    }
    MATRIX diagonal_deg_matrix = diagonalDegreeMatrix(weighted_adj_matrix, no_points);
    if (!diagonal_deg_matrix) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);
        freeMatrix(weighted_adj_matrix);
        return NULL;
    }
    MATRIX graph_laplacian = graphLaplacian(weighted_adj_matrix, diagonal_deg_matrix, no_points);
    if (!graph_laplacian) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        freeMatrix(c_point_list);
        freeMatrix(weighted_adj_matrix);
        freeMatrix(diagonal_deg_matrix);
        return NULL;
    }
    PyObject* py_gl = PyList_New(no_points);
    for (i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(py_gl, i, row);        //raises exception
        for (j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(graph_laplacian[i][j]));
        }
    }
    freeMatrix(c_point_list);
    freeMatrix(weighted_adj_matrix);
    freeMatrix(diagonal_deg_matrix);
    freeMatrix(graph_laplacian);
    return py_gl;
}
static PyObject* jacob(PyObject* self, PyObject* args) {
    PyObject* py_point_list;
    char* how;
    int requested_k;
    if (!PyArg_ParseTuple(args, "Os*n", &py_point_list,&how,&requested_k)) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    if (requested_k >= 0 && (strcmp(how, "sorted") != 0)) {
        PyErr_SetString(PyExc_Exception, "how should be \"sorted\" if k>=0");
        return NULL;
    }
    size_t no_points = PyObject_Length(py_point_list);
    PyObject* point = PyList_GetItem(py_point_list, 0);
    size_t dimension = PyObject_Length(point);
    if (no_points != dimension) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    size_t i, j;
    //allocate c matrix
    MATRIX c_point_list = allocateNonSquareMatrix(no_points, dimension);
    if (!c_point_list) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < dimension; i++) {
        point = PyList_GetItem(py_point_list, i);
        for (j = 0; j < dimension; j++) {
            c_point_list[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    EIGEN_VALUES_VECTORS* jacobis_res = jacobis(c_point_list, dimension, 100, 0.00001);
    if (!jacobis_res) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    size_t k = -1;
    if (strcmp(how,"sorted")==0) {
        sortEigenVectors(jacobis_res, dimension);
        k = find_k(jacobis_res->eigen_values, dimension);
    }
    int col_limit = dimension;
    if (requested_k > 0) {
        col_limit = requested_k;
    }
    else if (requested_k == 0)
        col_limit = k;

    PyObject* py_eigen_vectors = PyList_New(dimension);
    MATRIX c_eigen_vectors = jacobis_res->eigen_vectors_as_columns;
    for (i = 0; i < dimension; i++) {
        PyObject* row = PyList_New(col_limit);
        PyList_SetItem(py_eigen_vectors, i, row);        //raises exception
        for (j = 0; j < col_limit; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(c_eigen_vectors[i][j]));
        }
    }
    PyObject* py_eigen_values = PyList_New(dimension);
    VECTOR c_eigen_values = jacobis_res->eigen_values;
    for (i = 0; i < dimension; i++) {
        PyList_SetItem(py_eigen_values, i, PyFloat_FromDouble(c_eigen_values[i]));
    }

    freeMatrix(c_point_list);
    freeMatrix(jacobis_res->eigen_vectors);
    freeMatrix(jacobis_res->eigen_vectors_as_columns);
    free(jacobis_res->eigen_values);
    free(jacobis_res);
    return Py_BuildValue("OOn",py_eigen_values,py_eigen_vectors,k);
}


static PyMethodDef methods[] = {
    {"spk",(PyCFunction)fit, METH_VARARGS, PyDoc_STR("takes 2 python lists, max iteration value, Convergence value")},
    {"wam",(PyCFunction)wam,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the weight adjacency matrix")},
    {"ddg",(PyCFunction)ddg,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the diagonal degree matrix")},
    {"gl",(PyCFunction)gl,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the graph laplacian matrix")},
    {"jacobi",(PyCFunction)jacob,METH_VARARGS,PyDoc_STR("jacobis on a symmetric matrix,second argument should be \"sorted\" for spk purposes\n, returns(values,vectors matrix,k)")},
      /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     

};

static struct PyModuleDef mykmeanssp = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject* m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}