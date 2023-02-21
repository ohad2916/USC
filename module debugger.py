import mykmeanssp
matrix = [[2.1231, 123.123], [123.12334, 23.211], [264.5, 67778.4]]
symmetric_matrix = [[1, 0, 0], [0, 4, -1], [0, -1, 2]]
wam = mykmeanssp.wam(matrix)
ddg = mykmeanssp.ddg(matrix)
gl = mykmeanssp.gl(matrix)
values, vectors, k = mykmeanssp.jacobi(symmetric_matrix, "sorted", 0)
print(k)
print("")
print(','.join(f'{v:.4f}' for v in values))
print("")
for line in vectors:
    print(','.join(f'{c:.4f}' for c in line))
# for line in wam:
#     print(','.join(f'{c:.4f}' for c in line))
# print("")
# for line in ddg:
#     print(','.join(f'{c:.4f}' for c in line))
# print("")
# for line in gl:
#     print(','.join(f'{c:.4f}' for c in line))
# print("")
