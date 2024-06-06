import matplotlib.pyplot as plt

# 读取数据
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

# 从文件中读取数据
x1, y1 = read_data('fft_f1(1).txt')
x2, y2 = read_data('fft_f1(2).txt')

# 绘制图像
plt.figure(figsize=(12, 6))

plt.plot(x1, y1, marker='o', color='blue', label='n=16')
plt.plot(x2, y2, marker='o', color='red', label='n=128')

plt.title('FFT of f1')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

plt.tight_layout()
plt.show()