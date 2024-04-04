import numpy as np

# Define the A1 matrix
A1 = np.array([[1.0 / (9 - i - j) for j in range(5)] for i in range(5)])
print(A1)
# Compute the eigenvalues of A1
eigenvalues_A1 = np.linalg.eigvals(A1)

print("Eigenvalues of A1:", eigenvalues_A1)

# Define the A2 matrix
A2 = np.array([
    [4, -1, 1, 3],
    [16, -2, -2, 5],
    [16, -3, -1, 7],
    [6, -4, 2, 9]
])

# Compute the eigenvalues of A2
eigenvalues_A2 = np.linalg.eigvals(A2)

print("Eigenvalues of A2:", eigenvalues_A2)