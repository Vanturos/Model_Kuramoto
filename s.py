import numpy as np
import matplotlib.pyplot as plt

K = int(input('Введите коэффициент связи:'))

# Задание параметров модели
n = 4 # количество гребцов в экипаже
f0 = 1 # начальная частота опускания весла
sigma = 0.15 # стандартное отклонение для начальной частоты опускания весла
tmax = 50 # длительность симуляции
dt = 0.01 # шаг по времени
K1 = 1 # коэффициент связи между первым гребцом и остальными
K2 = 1/K # коэффициент связи между остальными гребцами

# Создание массивов для хранения частот опускания весла каждого гребца
f = np.zeros((n, int(tmax/dt)))

# Инициализация начальных частот опускания весла
f[:,0] = np.random.normal(f0, sigma*f0, n)

# Создание массивов для хранения частот опускания весла каждого гребца
f_2 = np.zeros((n, int(tmax / dt)))

# Инициализация начальных частот опускания весла
f_2[:, 0] = np.random.normal(f0, sigma * f0, n)

# Функция, определяющая изменение частот опускания весла для каждого гребца
def kuramoto1(f, K1, K2):
    dfdt = np.zeros(n)
    for i in range(n):
        for j in range(n):
            dfdt[i] += K1*(j>i)*np.sin(f[j]-f[i]) + K2*(j<i)*np.sin(f[j]-f[i])
    return dfdt

# Решение уравнения Курамото методом Эйлера
for i in range(len(f[0])-1):
    f[:,i+1] = f[:,i] + kuramoto1(f[:,i], K1, K2)*dt

print('Для случая все - все:')
# Поиск критической точки
for i in range(len(f[0]) - 1):
    a, b, c, d = f[:, i]
    x = 0.001
    if abs(a-b)<x and abs(b-c)<x and abs(c-d)<x and abs(d-a)<x:
        print("Критическая точка = {}\nЧастота в ней = {}".format(round(dt*(i+1),2), round(a, 2)))
        break

# Функция, определяющая изменение частот опускания весла для каждого гребца
def kuramoto2(f, K1):
    dfdt = np.zeros(n)
    for i in range(n):
        for j in range(n):
            dfdt[i] += K1*(j==0)*np.sin(f[j]-f[i])
    return dfdt

# Решение уравнения Курамото методом Эйлера
for i in range(len(f_2[0]) - 1):
    f_2[:, i + 1] = f_2[:, i] + kuramoto2(f_2[:, i], K1) * dt

print('\nДля случая один - все:')
# Поиск критической точки
for i in range(len(f_2[0]) - 1):
    a, b, c, d = f_2[:, i]
    x = 0.001
    if abs(a-b)<x and abs(b-c)<x and abs(c-d)<x and abs(d-a)<x:
        print("Критическая точка = {}\nЧастота в ней = {}".format(round(dt*(i+1),2), round(a, 2)))
        break

# Визуализация результатов
fig,(ax1, ax2) = plt.subplots(2)
ax1.set_title("График зависимости частоты от времени для случая все - все")
ax2.set_title("График зависимости частоты от времени для случая один - все")
fig.set_size_inches(18.5, 10.5)
plt.subplots_adjust(wspace=0.3, hspace=0.3)
for i in range(n):
    ax1.plot(np.arange(0,tmax,dt), f[i,:], label='гребец {}'.format(i+1))
for i in range(n):
    ax2.plot(np.arange(0,tmax,dt), f_2[i, :], label='гребец {}'.format(i + 1))
ax1.set(xlabel='Время, с', ylabel='Частота опускания весла, 1/t')
ax2.set(xlabel='Время, с', ylabel='Частота опускания весла, 1/t')
plt.legend()
plt.show()


