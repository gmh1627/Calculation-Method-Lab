import numpy as np
from scipy.integrate import quad
from math import sin, sqrt, log
import matplotlib.pyplot as plt

# 定义加速度函数
def ax(t):
    return sin(t) / (1 + sqrt(t))

def ay(t):
    return log(t + 1) / (t + 1)

# 初始化速度和位移数组
vx = np.zeros(101)
vy = np.zeros(101)
x = np.zeros(101)
y = np.zeros(101)

# 计算每个时间点的速度和位移
for i in range(1, 101):
    t = i / 10.0
    vx[i], _ = quad(ax, 0, t)
    vy[i], _ = quad(ay, 0, t)
    print(f"At t={t:.1f}s, vx={vx[i]:.3f} m/s, vy={vy[i]:.3f} m/s")

# 计算位移
for i in range(1, 101):
    t = i / 10.0
    x[i], _ = quad(lambda t: vx[int(t*10)], 0, t)
    y[i], _ = quad(lambda t: vy[int(t*10)], 0, t)
    print(f"At t={t:.1f}s, x={x[i]:.3f} m, y={y[i]:.3f} m")

# 创建时间数组
time = np.linspace(0, 10, 101)

# 绘制速度图
plt.figure(figsize=(12, 6))
plt.subplot(2, 2, 1)
plt.plot(time, vx, label='vx')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Velocity in x direction')
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(time, vy, label='vy')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Velocity in y direction')
plt.grid(True)

# 绘制位移图
plt.subplot(2, 2, 3)
plt.plot(time, x, label='x')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Displacement in x direction')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(time, y, label='y')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Displacement in y direction')
plt.grid(True)

# 显示图像
plt.tight_layout()
plt.show()