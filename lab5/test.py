import numpy as np

# 定义矩阵A和向量b
A = np.array([[1, 4, -3], [1, 4, -3], [1, 4, -3]])
b = np.array([1, 1, 2])

# 使用numpy.linalg.solve解决Ax = b
x = np.linalg.solve(A, b)

# 打印结果
print(x)