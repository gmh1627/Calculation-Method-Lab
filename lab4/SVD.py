import numpy as np
from scipy.linalg import svd

# 定义矩阵
A = np.array([[0.726249, 0.915339, 0.594324],
              [0.515358, 0.975149, 0.661561],
              [0.528652, 0.788493, 0.0741007],
              [0.32985, 0.583744, 0.0309722]])

# 进行奇异值分解
U, s, VT = svd(A)

print("U:\n", U)
print("s:\n", s)
print("VT:\n", VT)