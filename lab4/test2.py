import numpy as np

A2 = np.array([
    [1.71851, 1.66005, 1.14971, 0.792285],
    [1.66005, 1.65417, 1.09037, 0.759718],
    [1.14971, 1.09037, 0.906686, 0.636949],
    [0.792285, 0.759718, 0.636949, 0.450517],
])

# 创建一个全为0的向量b
b = np.zeros(A2.shape[0])

# 创建一个A的副本
A = np.copy(A2)

# 从A的对角线元素中减去4.52037
np.fill_diagonal(A, np.diagonal(A) - 4.52037)

print(A)
print(np.linalg.det(A))
# 解线性方程Ax=b
x = np.linalg.solve(A, b)

print(x)