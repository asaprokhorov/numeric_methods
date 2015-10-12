import matplotlib.pyplot as plt
import numpy
import time
import math


def f(x):
    return numpy.tan(x / 2)


a = -1.0
b = 2.0
N = 20
h = (b - a) / N

xFunc = []
yFunc = []
differences = []

for i in range(N + 1):
    xFunc.append(a + i * h)
    yFunc.append(f(a + i * h))


def finite_diff2(powx, i):
    diff = 0.0
    for j in range(powx + 1):
        diff += ((-1)**(powx - j))*(math.factorial(powx) / (math.factorial(j) * math.factorial((powx - j)) * f(xFunc[i + j])))
    return diff


def newton_forward2(x):
    result = f(a)
    elem = 1
    t = (x - a) / h
    for i in range(1, N + 1):
        elem *= (t - i + 1) / i
        result += elem * finite_diff2(i, 0)
    return result


xNewt2 = []
yNewt2 = []
t = time.clock()
for i in range(N + 1):
    xNewt2.append(a + i * h)
    yNewt2.append(newton_forward2(a + i * h))
print("Time spent on newton(Differences calculated)")
print(time.clock() - t)


def newton_forward(x, differences2):
    result = f(a)
    elem = 1
    t = (x - a) / h
    for i in range(1, N + 1):
        elem *= (t - i + 1) / i
        result += elem * differences2[i][0]
    return result


xNewt = []
yNewt = []

t = time.clock()
differences.append(yFunc)
for i in range(1, N + 1):
    differences.append([])
    for j in range(0, N - i + 1):
        differences[i].append(differences[i - 1][j + 1] - differences[i - 1][j])

for i in range(N + 1):
    xNewt.append(a + i * h)
    yNewt.append(newton_forward(a + i * h, differences))
print("Time spent on newton(Differences from array)")
print(time.clock() - t)


n = N // 2


def gauss_backward(x):
    x0 = xFunc[n]
    t = (x - x0) / h
    result = f(x0)
    elem = 1
    for i in range(1, n + 1):
        elem *= (t - i + 1) / (2 * i - 1)
        result += elem * differences[2 * i - 1][n - i]
        elem *= (t + i) / (2 * i)
        result += elem * differences[2 * i][n - i]
    elem *= (t - n) / (N)
    result += elem * differences[N][0]
    return result


xGauss = []
yGauss = []

for i in range(N + 1):
    xGauss.append(a + i * h)
    yGauss.append(gauss_backward(a + i * h))

print("Newton difference")
for i in range(N + 1):
    print(i, abs(yFunc[i] - yNewt[i]))

print("Gauss difference")
for i in range(N + 1):
    print(i, abs(yFunc[i] - yGauss[i]))

# Graphic f(x)
plt.plot(xFunc, yFunc)

# Newton Graphic
plt.plot(xNewt, yNewt)

# Gauss Graphic
plt.plot(xGauss, yGauss)

plt.show()