import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

def readIrisData(filename):
    data = np.genfromtxt(filename, delimiter=',', dtype='float', encoding=None)
    return data[:, :4].T, data[:, 4]

X, labels = readIrisData("iris.txt")
n1, n2 = X.shape
X -= X.mean(axis=1, keepdims=True)
XT = X.T
XXT = np.dot(X, XT)

Var = XXT / n2
x, EigenVector = eigh(Var)
x = sorted(x, reverse=True)

P = EigenVector[:, -2:].T
Y = np.dot(P, X)

plt.figure()
label_set = set(labels)
colors = ['r', 'g', 'b']
shapes = ['o', 's', '^']

for i, label in enumerate(label_set):#enumerate函数返回每个标签及其索引
    x = [Y[0, j] for j in range(Y.shape[1]) if labels[j] == label]
    y = [Y[1, j] for j in range(Y.shape[1]) if labels[j] == label]
    plt.scatter(x, y, color=colors[i], marker=shapes[i], label=label)

plt.legend()#添加图例
plt.show()