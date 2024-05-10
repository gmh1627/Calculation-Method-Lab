import matplotlib.pyplot as plt

with open('trajectory.txt', 'r') as file:
    lines = file.readlines()

x, y = zip(*[(float(line.split()[0]), float(line.split()[1])) for line in lines])

plt.plot(x, y, 'o-')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Trajectory of the particle')
plt.show()