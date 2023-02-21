import mykmeanssp as mk
import pandas as pd
import numpy as np
import sys
import math

arg_c = len(sys.argv)
np.random.seed(0)
iter_ = 200
epsilon = 0.001
k = -1

if arg_c > 4 or arg_c < 3:
    print("An Error Has Occurred")
    sys.exit()
try:
    # total arguments
    if arg_c == 3:
        goal = sys.argv[1]
        file_name = sys.argv[2]
    elif arg_c == 4:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]

except Exception as e:
    print("An Error Has Occurred")
    sys.exit()

if arg_c > 4 or arg_c < 3:
    print("An Error Has Occurred")
    sys.exit()


def euc_d(p1, p2, dim):
    np_arr_p1 = np.array(p1)
    np_arr_p2 = np.array(p2)

    sum_ = 0.0
    for i in range(dim):
        sum_ += (np_arr_p1[i] - np_arr_p2[i]) ** 2

    return math.sqrt(sum_)


def init_centroids(data, num_clusters):
    centroids = []
    centroids_indices = []
    prob = np.zeros(data.shape[0])

    first_centroid_index = int(np.random.choice(data.shape[0]))
    first_centroid = data.iloc[first_centroid_index]
    centroids.append(first_centroid)
    centroids_indices.append(first_centroid_index)

    for i in range(num_clusters - 1):  # init num_clusters-1 more centroids
        sum_dx = 0
        for j in range(data.shape[0]):  # calculate dx for each point
            p = data.iloc[j]
            # we first calculate the distance to the last centroid and then check if there is a closer one.
            dx = euc_d(centroids[len(centroids) - 1], p, len(p) - 1)
            for m in range(len(centroids) - 1):
                dist_to_next_centroid = euc_d(centroids[m], p, len(p) - 1)
                dx = min(dx, dist_to_next_centroid)

            prob[j] = dx
            sum_dx += dx

        prob /= sum_dx

        next_centroid_index = int(np.random.choice(data.shape[0], p=prob))
        new_centroid = data.iloc[next_centroid_index]
        centroids.append(new_centroid)
        centroids_indices.append(next_centroid_index)

    return centroids, centroids_indices


if goal == 'spk':
    X_df = pd.read_csv(file_name, header=None)
    X = X_df.values.tolist()
    L = mk.gl(X)  # graph Laplacian

    if k == -1:
        reg_jacobi_values, U, reg_k = mk.jacobi(L, "sorted", 0)
        starting_centroids, starting_centroids_indices = init_centroids(pd.DataFrame(U), reg_k)

    elif k >= 0:
        reg_jacobi_values, U, reg_k = mk.jacobi(L, "sorted", k)
        starting_centroids, starting_centroids_indices = init_centroids(pd.DataFrame(U), k)

    try:
        kmeans_res = mk.spk(X, starting_centroids, iter_, epsilon)
    except:
        print("An Error Has Occurred")
        sys.exit()

    print(",".join(f'{i:.0f}' for i in starting_centroids_indices))
    for cluster in kmeans_res:
        print(','.join(f'{c:.4f}' for c in cluster))


if goal in ['wam', 'ddg', 'gl']:
    X_df = pd.read_csv(file_name, header=None)
    X = X_df.to_list()

    if goal == 'wam':
        Matrix = mk.wam(X)  # weighted adjacency matrix
    if goal == 'ddg':
        Matrix = mk.ddg(X)  # diagonal degree matrix
    if goal == 'gl':
        Matrix = mk.gl(X)  # graph Laplacian

    for line in Matrix:
        print(",".join(f'{val:.4f}' for val in line))

if goal == 'jacobi':
    sym_Matrix_df = pd.read_csv(file_name, header=None)
    sym_Matrix = sym_Matrix_df.to_list()

    reg_jacobi_values, reg_jacobi_vectors, reg_k = mk.jacobi(sym_Matrix, "unsorted", -1)

    print(','.join(f'{c:.4f}' for c in reg_jacobi_values))
    for line in reg_jacobi_vectors:
        print(','.join(f'{val:.4f}' for val in line))
