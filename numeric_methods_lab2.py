import matplotlib.pyplot as plt
import numpy


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

differences.append(yFunc)
for i in range(1, N + 1):
    differences.append([])
    for j in range(0, N - i + 1):
        differences[i].append(differences[i - 1][j + 1] - differences[i - 1][j])


def newton_forward(x):
    result = f(a)
    elem = 1
    t = (x - a) / h
    for i in range(1, N + 1):
        elem *= (t - i + 1) / i
        result += elem * differences[i][0]
    return  result


xNewt = []
yNewt = []

for i in range(N + 1):
    xNewt.append(a + i * h)
    yNewt.append(newton_forward(a + i * h))


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
    return result


xGauss = []
yGauss = []

for i in range(N + 1):
    xGauss.append(a + i * h)
    yGauss.append(gauss_backward(a + i * h))

plt.plot(xFunc, yFunc)
plt.plot(xNewt, yNewt)
plt.plot(xGauss, yGauss)

plt.show()