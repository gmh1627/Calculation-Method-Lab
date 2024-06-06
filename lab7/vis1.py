import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# 定义加速度函数
def ax(t):
    return np.sin(t) / (1 + np.sqrt(t))

def ay(t):
    return np.log(t + 1) / (t + 1)

# 定义时间间隔
t = np.arange(0, 10.1, 0.1)

# 数值积分以计算速度
vx = integrate.cumtrapz(ax(t), t, initial=0)
vy = integrate.cumtrapz(ay(t), t, initial=0)

# 数值积分以计算位移
sx = integrate.cumtrapz(vx, t, initial=0)
sy = integrate.cumtrapz(vy, t, initial=0)

# 打印结果
print("时间点 t:")
print(t)
print("x 方向的速度 vx:")
print(vx)
print("x 方向的位移 sx:")
print(sx)
print("y 方向的速度 vy:")
print(vy)
print("y 方向的位移 sy:")
print(sy)

# 可视化结果
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t, vx, label='vx (x方向速度)')
plt.xlabel('时间 (s)')
plt.ylabel('速度 (m/s)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(t, sx, label='sx (x方向位移)')
plt.xlabel('时间 (s)')
plt.ylabel('位移 (m)')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, vy, label='vy (y方向速度)')
plt.xlabel('时间 (s)')
plt.ylabel('速度 (m/s)')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(t, sy, label='sy (y方向位移)')
plt.xlabel('时间 (s)')
plt.ylabel('位移 (m)')
plt.legend()

plt.tight_layout()
plt.show()
