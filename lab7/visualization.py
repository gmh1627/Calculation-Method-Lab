import matplotlib.pyplot as plt

with open('trajectory.txt', 'r') as file:
    lines = file.readlines()

x, y = zip(*[(float(line.split()[0]), float(line.split()[1])) for line in lines])

plt.plot(x, y, 'o-')
plt.xlabel('X')  # Change the x-axis label to 'X'
plt.ylabel('Y')  # Change the y-axis label to 'Y'
plt.title('Plot of points with smooth curve')
plt.show()