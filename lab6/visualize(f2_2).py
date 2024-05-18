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
x_discrete, y_discrete = read_data('discrete_f2.txt')
x_ifft, y_ifft = read_data('ifft_f2(2).txt')

# 绘制图像
plt.figure(figsize=(12, 6))

plt.plot(x_discrete, y_discrete, marker='o', color='blue', label='Discrete f1')
plt.plot(x_ifft, y_ifft, marker='o', color='red', label='IFFT of f1')

plt.title('Discrete and IFFT of f1')
plt.xlabel('x')
plt.ylabel('f1')
plt.legend()

plt.tight_layout()
plt.show()