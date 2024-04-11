import numpy as np

# Define the A2 matrix

A2 = np.array([
    [1.71851, 1.66005, 1.14971, 0.792285],
    [1.66005, 1.65417, 1.09037, 0.759718],
    [1.14971, 1.09037, 0.906686, 0.636949],
    [0.792285, 0.759718, 0.636949, 0.450517],
])
"""
A2 = np.array([
    [1, 2, 1],
    [1, 2, 2],
    [3, 8, 3]
])
"""
# Compute the eigenvalues of A2
eigenvalues_A2, eigenvectors_A2 = np.linalg.eig(A2)

print("Eigenvalues of A2:", eigenvalues_A2)
print("Eigenvectors of A2:")
for i in range(len(eigenvalues_A2)):
    print("Eigenvector for eigenvalue {}: {}".format(eigenvalues_A2[i], eigenvectors_A2[:, i]))