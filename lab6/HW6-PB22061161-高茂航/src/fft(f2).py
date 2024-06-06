import matplotlib.pyplot as plt

def read_data(file_path):
    with open(file_path, 'r') as file:
        data = file.readlines()
    x = []
    y = []
    for line in data:
        values = line.split()
        x.append(float(values[0]))
        y.append(float(values[1]))
    return x, y

x, y = read_data('fft_f2.txt')

# 绘制图像
plt.plot(x, y, marker='o', label='FFT of f2')

plt.title('FFT of f2')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

plt.show()