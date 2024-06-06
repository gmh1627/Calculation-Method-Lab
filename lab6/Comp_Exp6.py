# import numpy as np
# import matplotlib.pyplot as plt

# # 定义函数f(t)
# def f(t):
#     return 0.7 * np.sin(2 * np.pi * 2 * t) + np.sin(2 * np.pi * 5 * t)

# # 采样时间间隔
# dt = 1 / 16

# # 采样时间范围
# t_values = np.arange(0, 1, dt)

# # 计算向量f
# f_values = f(t_values)
# print(f_values)

# # 执行快速傅立叶变换
# g_values = np.fft.fft(f_values)
# print(g_values)

# # 取一半的向量g
# half_g_values = g_values[:len(g_values)//2]

# # 计算每个向量的模长并乘以2
# magnitudes = np.abs(half_g_values) * 2

# # 绘制频率域上能量分布图
# plt.plot(np.arange(0, 0.5, dt), magnitudes)
# plt.xlabel('Frequency')
# plt.ylabel('Magnitude')
# plt.title('Frequency Domain Energy Distribution')
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

# # 从文件中读取复向量数据
# def read_complex_vector(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
#     complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
#     return np.array(complex_vector)

# # 执行逆傅里叶变换
# def perform_ifft(vector):
#     return np.fft.ifft(np.fft.fft(vector))

# # 绘制图像
# def plot_vectors(original, transformed, n, subplot_position):
#     t_values = np.array([(i-1)/n for i in range(1, n+1)])
#     plt.subplot(2, 1, subplot_position)
#     plt.plot(t_values, original.real, label='Original Vector')
#     plt.plot(t_values, transformed.real, label='IFFT(FFT(Vector))', linestyle='--')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Amplitude')
#     plt.legend()
#     plt.title(f'n = {n}')

# # 主函数
# def main():
#     filenames = ['C:\\Users\\23103\\Desktop\\f1_1.txt', 'C:\\Users\\23103\\Desktop\\h1_1.txt']  # 替换为你的文件名
#     n_values = [16, 16]  # 对应文件中的n值

#     plt.figure(figsize=(12, 8))
    
#     for i, (filename, n) in enumerate(zip(filenames, n_values), 1):
#         original_vector = read_complex_vector(filename)
#         transformed_vector = perform_ifft(original_vector)
#         plot_vectors(original_vector, transformed_vector, n, i)
    
#     plt.tight_layout()
#     plt.show()

# if __name__ == '__main__':
#     main()

# import numpy as np
# import matplotlib.pyplot as plt

# # 从文件中读取复向量数据
# def read_complex_vector(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
#     complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
#     return np.array(complex_vector)

# # 执行逆傅里叶变换
# def perform_ifft(vector):
#     return np.fft.ifft(np.fft.fft(vector))

# # 绘制图像
# def plot_vectors(original1, transformed1, n1, original2, transformed2, n2):
#     t_values1 = np.array([(i-1)/n1 for i in range(1, n1+1)])
    
#     plt.figure(figsize=(12, 8))
    
#     plt.plot(t_values1, original1.real, label=f'Original Vector n={n1}')
#     plt.plot(t_values1, transformed1.real, label=f'IFFT(FFT(Vector)) n={n2}', linestyle='--')
    
#     plt.xlabel('Time (s)')
#     plt.ylabel('Amplitude')
#     plt.legend()
#     plt.title('Original and IFFT(FFT) Vectors')
#     plt.show()

# # 主函数
# def main():
#     filenames = ['C:\\Users\\23103\\Desktop\\f1_1.txt', 'C:\\Users\\23103\\Desktop\\h1_1.txt']  # 替换为你的文件名
#     n_values = [16, 16]  # 对应文件中的n值

#     original_vector1 = read_complex_vector(filenames[0])
#     transformed_vector1 = read_complex_vector(filenames[1])
    
#     plot_vectors(original_vector1, transformed_vector1, n_values[0], n_values[1])

# if __name__ == '__main__':
#     main()

# import numpy as np
# import matplotlib.pyplot as plt

# # 从文件中读取复向量数据
# def read_complex_vector(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
#     complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
#     return np.array(complex_vector)

# # 绘制图像
# def plot_vectors(original, transformed, n):
#     t_values = np.array([(i-1)/n for i in range(1, n+1)])
    
#     plt.figure(figsize=(12, 8))
    
#     plt.plot(t_values, original.real, label='Original Vector')
#     plt.plot(t_values, transformed.real, label='IFFT(FFT(Vector))', linestyle='--')
    
#     plt.xlabel('Time (s)')
#     plt.ylabel('Amplitude')
#     plt.legend()
#     plt.title(f'Original and IFFT(FFT) Vectors for n={n}')
#     plt.show()

# # 主函数
# def main():
#     original_filename = 'C:\\Users\\23103\\Desktop\\f2.txt'  # 替换为你的原向量文件名
#     transformed_filename = 'C:\\Users\\23103\\Desktop\\h2_1.txt'  # 替换为你的逆傅里叶变换结果文件名
#     transformed_filename2 = 'C:\\Users\\23103\\Desktop\\h2_2.txt'  # 替换为你的逆傅里叶变换结果文件名
#     n = 128  # 替换为你的n值

#     original_vector = read_complex_vector(original_filename)
#     transformed_vector = read_complex_vector(transformed_filename)
#     transformed_vector2 = read_complex_vector(transformed_filename2)
    
#     plot_vectors(original_vector, transformed_vector, transformed_vector2, n)

# if __name__ == '__main__':
#     main()

# import numpy as np
# import matplotlib.pyplot as plt

# # 从文件中读取复向量数据
# def read_complex_vector(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
#     complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
#     return np.array(complex_vector)

# # 绘制图像
# def plot_vectors(original, transformed, additional, n):
#     t_values = np.array([(i-1)/n for i in range(1, n+1)])
    
#     plt.figure(figsize=(12, 8))
    
#     plt.plot(t_values, original.real, label='Original Vector')
#     plt.plot(t_values, transformed.real, label='IFFT(FFT(Vector))', linestyle='--')
#     plt.plot(t_values, additional.real, label='Additional Result', linestyle=':')
    
#     plt.xlabel('Time (s)')
#     plt.ylabel('Amplitude')
#     plt.legend()
#     plt.title(f'Original, IFFT(FFT) Vectors and Additional Result for n={n}')
#     plt.show()

# # 主函数
# def main():
#     original_filename = 'C:\\Users\\23103\\Desktop\\f2.txt'  # 替换为你的原向量文件名
#     transformed_filename = 'C:\\Users\\23103\\Desktop\\h2_1.txt'  # 替换为你的逆傅里叶变换结果文件名
#     additional_filename = 'C:\\Users\\23103\\Desktop\\h2_2.txt'  # 替换为你的额外结果文件名
#     n = 128  # 替换为你的n值

#     original_vector = read_complex_vector(original_filename)
#     transformed_vector = read_complex_vector(transformed_filename)
#     additional_vector = read_complex_vector(additional_filename)
    
#     plot_vectors(original_vector, transformed_vector, additional_vector, n)

# if __name__ == '__main__':
#     main()

# import numpy as np
# import matplotlib.pyplot as plt

# # 从文件中读取复向量数据
# def read_complex_vector(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()
#     complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
#     return np.array(complex_vector)

# # 绘制频谱图
# def plot_frequency_spectrum(frequencies, magnitudes, n):
#     plt.figure(figsize=(10, 6))
#     plt.plot(frequencies, magnitudes, label=f'n={n}')
#     plt.xlabel('Frequency (Hz)')
#     plt.ylabel('|G(f)|')
#     plt.title(f'Frequency Spectrum for n={n}')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# # 主函数
# def main():
#     filename_n16 = 'C:\\Users\\23103\\Desktop\\g1_1.txt'  # 替换为n=16的文件名
#     filename_n128 = 'C:\\Users\\23103\\Desktop\\g1_2.txt'  # 替换为n=128的文件名

#     # 读取复向量数据
#     g_n16 = read_complex_vector(filename_n16)
#     g_n128 = read_complex_vector(filename_n128)

#     # 计算幅值
#     magnitudes_n16 = np.abs(g_n16)
#     magnitudes_n128 = np.abs(g_n128)

#     # 计算频率轴
#     n16 = 16
#     n128 = 128
#     sample_rate = n16  # 采样率假设为样本数，因为时域范围为0~1
#     freqs_n16 = np.fft.fftfreq(n16, d=1/sample_rate)
#     freqs_n128 = np.fft.fftfreq(n128, d=1/sample_rate)

#     # 只取前半部分（因为频谱是对称的）
#     half_n16 = n16 // 2
#     half_n128 = n128 // 2
#     freqs_n16 = freqs_n16[:half_n16]
#     magnitudes_n16 = magnitudes_n16[:half_n16]
#     freqs_n128 = freqs_n128[:half_n128]
#     magnitudes_n128 = magnitudes_n128[:half_n128]

#     # 绘制频谱图
#     plt.figure(figsize=(12, 8))
#     plt.plot(freqs_n16, magnitudes_n16, label=f'n=16', marker='o')
#     plt.plot(freqs_n128, magnitudes_n128, label=f'n=128', marker='x')
#     plt.xlabel('Frequency (Hz)')
#     plt.ylabel('|G(f)|')
#     plt.title('Frequency Spectrum for different n')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# if __name__ == '__main__':
#     main()

import numpy as np
import matplotlib.pyplot as plt

# 从文件中读取复向量数据
def read_complex_vector(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    complex_vector = [complex(*map(float, line.strip()[1:-1].split(','))) for line in lines]
    return np.array(complex_vector)

# 绘制频谱图
def plot_frequency_spectrum(frequencies, magnitudes, n):
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, magnitudes, label=f'n={n}', marker='o')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('|G(f)|')
    plt.title(f'Frequency Spectrum for n={n}')
    plt.legend()
    plt.grid(True)
    plt.show()

# 主函数
def main():
    filename_n16 = 'C:\\Users\\23103\\Desktop\\g2.txt'  # 替换为n=16的文件名

    # 读取复向量数据
    g_n16 = read_complex_vector(filename_n16)

    # 计算幅值
    magnitudes_n16 = np.abs(g_n16)

    # 计算频率轴
    n16 = 16
    sample_rate = n16  # 采样率假设为样本数，因为时域范围为0~1
    freqs_n16 = np.fft.fftfreq(n16, d=1/sample_rate)

    # 只取前半部分（因为频谱是对称的）
    half_n16 = n16 // 2
    freqs_n16 = freqs_n16[:half_n16]
    magnitudes_n16 = magnitudes_n16[:half_n16]

    # 绘制频谱图
    plot_frequency_spectrum(freqs_n16, magnitudes_n16, n16)

if __name__ == '__main__':
    main()
