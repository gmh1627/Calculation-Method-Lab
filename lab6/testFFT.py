import numpy as np
import matplotlib.pyplot as plt

# 定义函数
def f1(t):
    return 0.7 * np.sin(2 * np.pi * 2 * t) + np.sin(2 * np.pi * 5 * t)

# 设置采样率和时间范围
sample_rates = [128]  # 只使用抽样频率为128的情况
time_duration = 1  # 总时间范围

# 对每个采样率进行FFT并绘图
for sample_rate in sample_rates:
    # 生成采样点
    time_range = np.arange(0, time_duration, 1/sample_rate)
    samples = f1(time_range)

    # 计算FFT
    fft_result = np.fft.fft(samples)
    fft_result = np.fft.fftshift(fft_result)  # Shift zero freq to center

    # 计算频率
    freqs = np.fft.fftfreq(len(samples), d=1/sample_rate)
    freqs = np.fft.fftshift(freqs)  # Shift zero freq to center

    # 计算IFFT
    ifft_result = np.fft.ifft(np.fft.ifftshift(fft_result))

    # 创建一个新的图形窗口
    plt.figure(figsize=(10, 5))

    # 在同一图上绘制f1的离散点折线图和IFFT图像
    plt.plot(time_range, samples, label='f1 discrete points')
    plt.plot(time_range, ifft_result.real, label='IFFT result')
    plt.title(f'f1 discrete points and IFFT result with sample rate {sample_rate}')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.grid(True) # 显示网格
    plt.legend()

    # 创建一个新的图形窗口
    plt.figure(figsize=(10, 5))

    # 绘制FFT图像
    plt.plot(freqs, np.abs(fft_result), label='FFT result')
    plt.title(f'FFT of f1 with sample rate {sample_rate}')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.legend()

    # 显示图形
    plt.show()