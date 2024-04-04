import numpy as np

# Define the A2 matrix
A2 = np.array([
    [1, -1, 0],
    [-1, 2, 2],
    [0, 2, 3],
])

# Compute the eigenvalues of A2
eigenvalues_A2 = np.linalg.eigvals(A2)

print("Eigenvalues of A2:", eigenvalues_A2)