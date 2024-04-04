import numpy as np

def solve_linear_system(epsilon):
    n = 100
    a = 1.0 / 2
    h = 1.0 / n

    # 初始化A和b
    A = np.zeros((n-1, n-1))
    b = np.full((n-1,), a * h * h)

    # 填充A
    for i in range(n - 2):
        A[i][i + 1] = epsilon + h
        A[i + 1][i] = epsilon
    for i in range(n - 1):
        A[i][i] = -(2 * epsilon + h)

    # 更新b
    b[98] -= epsilon + h

    # 解方程
    y = np.linalg.solve(A, b)

    return y

def precise_sol(x, epsilon):
    a = 1.0 / 2
    y = (1 - a) * (1 - np.exp(-x / epsilon)) / (1 - np.exp(-1 / epsilon)) + a * x
    return y

def calculate_error(y, presol):
    error = np.abs(y - presol).mean()
    return error

# 使用函数
for epsilon in [1, 0.1, 0.01, 0.0001]:
    y = solve_linear_system(epsilon)
    x = np.linspace(0, 1, 99)
    presol = precise_sol(x, epsilon)
    error = calculate_error(y, presol)
    print(f"epsilon={epsilon}时，相对误差为{error}")