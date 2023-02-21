import mykmeanssp as mk

# example matrix : 1 0 0
#                  0 4 -1
#                  0,-1,2
matrix = [[1, 0, 0], [0, 4, -1], [0, -1, 2]]

# wam - weighted adjacency matrix. example:
wam = mk.wam(matrix)

# ddg - diagonal degree matrix, does not need the wam. example:
ddg = mk.ddg(matrix)

# gl - computes the graph laplacian. example:
gl = mk.gl(matrix)

# jacobi - provides the eigen values,vectors(as columns!!) of a symmetric matrix.
#          first argument is the matrix (should be the graph laplacian in case spk)
#          second is a string (how) - if "sorted", the eigen values,vectors returned will be
#          in ascending order also third return value is the K computed by the eigen-gap
#          heuristic, only in case how == "sorted" . third argument is the k inputted by
#          user, default should be ( == 0) and it will return the K first eigenvectors
#          where K is computed by the eigen-gap heuristic. third argument should be >=0
#          only if how = "sorted". examples:
# computing regular jacobi (not case 'spk'):
reg_jacobi_values, reg_jacobi_vectors, reg_k = mk.jacobi(matrix, "unsorted", -1)

print(','.join(f'{c:.3f}' for c in reg_jacobi_values))
for line in reg_jacobi_vectors:
    print(','.join(f'{c:.4f}' for c in line))
print(f'computed K is: {reg_k} ,defaults in case "unsorted" \n')
# computing full jacobi (case 'spk') with the K from eigen-gap heuristic:
full_jacobi_values, full_jacobi_vectors, computed_K = mk.jacobi(matrix, "sorted", 0)
print(','.join(f'{c:.3f}' for c in full_jacobi_values))
for line in full_jacobi_vectors:
    print(','.join(f'{c:.4f}' for c in line))
print(f'computed K is: {computed_K}\n')
# same thing with K provided by user:
user_K = 2
user_jacobi_values, user_jacobi_vectors, reg_k2 = mk.jacobi(matrix, "sorted", user_K)
print(','.join(f'{c:.3f}' for c in user_jacobi_values))
for line in user_jacobi_vectors:
    print(','.join(f'{c:.4f}' for c in line))
print(f'computed K is: {reg_k2} ,but user asked for 2\n')
