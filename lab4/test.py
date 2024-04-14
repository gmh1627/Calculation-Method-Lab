import numpy as np
import math 

# Define the A2 matrix
A2 = np.array([
    [1.07251, 1.58415, 0.811742],
    [1.58415, 2.41048, 1.24756],
    [0.811742, 1.24756, 0.796375],
])

# Compute the eigenvalues of A2
eigenvalues_A2, eigenvectors_A2 = np.linalg.eig(A2)

print("sqrt of Eigenvalues of A2:", np.sqrt(eigenvalues_A2))
print("Eigenvectors of A2:")
for i in range(len(eigenvalues_A2)):
    print("Eigenvector for eigenvalue {}: {}".format(eigenvalues_A2[i], eigenvectors_A2[:, i]))