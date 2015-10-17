import matplotlib.pyplot as plt
import numpy


def phi(x, i):
    return x ** i


def f(x):
    return numpy.tan(x / 2)


# Input data
a = -1.0
b = 2.0
N = 20
n = N // 2
h = (b - a) / N
xFunc = []
yFunc = []
for i in range(N + 1):
    xFunc.append(a + i * h)
    yFunc.append(f(a + i * h))
differences = [yFunc]
for i in range(1, N + 1):
    differences.append([])
    for j in range(N - i + 1):
        differences[i].append(differences[i - 1][j + 1] - differences[i - 1][j])


# Newton forward method
def newton_forward(x):
    result = f(xFunc[0])
    elem = 1
    t = (x - xFunc[0]) / h
    for i in range(1, N + 1):
        elem *= (t - i + 1) / i
        result += elem * differences[i][0]
    return result


# Gauss backward method
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
    elem *= (t - n) / N
    result += elem * differences[N][0]
    return result


# Newton forward data
xNewt = []
yNewt = []
for i in range(N + 1):
    xNewt.append(a + i * h)
    yNewt.append(newton_forward(a + i * h))

# Gauss backward data
xGauss = []
yGauss = []
for i in range(N + 1):
    xGauss.append(a + i * h)
    yGauss.append(gauss_backward(a + i * h))
# Graphic f(x)
plt.plot(xFunc, yFunc)

# Newton Graphic
plt.plot(xNewt, yNewt)

# Gauss Graphic
plt.plot(xGauss, yGauss)


def scalar_phi(i, j):
    result = 0
    for k in range(N + 1):
        result += phi(xFunc[k], i) * phi(xFunc[k], j)
    return result


def scalar(i):
    result = 0
    for k in range(N + 1):
        result += f(xFunc[k]) * phi(xFunc[k], i)
    return result


matrix = []
for i in range(N + 1):
    matrix.append([])
    for j in range(N + 1):
        matrix[i].append(scalar_phi(i, j))
b = []
for i in range(N + 1):
    b.append(scalar(i))


class JordanGauss:
    def __init__(self, matr, col):
        self.matrix = matr
        self.b = col

    def make_t(self, k):
        t = []
        for i in range(N + 1):
            t.append([])
            for j in range(N + 1):
                if i == j:
                    t[i].append(1)
                else:
                    t[i].append(0)
        for i in range(N + 1):
            t[i][k] = -self.matrix[i][k] / self.matrix[k][k]
        t[k][k] = 1. / self.matrix[k][k]
        return t

    def multiply(self, matr):
        result = []
        for i in range(N + 1):
            result.append([])
            for j in range(N + 1):
                sum = 0
                for k in range(N + 1):
                    sum += self.matrix[k][j] * matr[i][k]
                result[i].append(sum)
        self.matrix = result
        result = []
        for i in range(N + 1):
            sum = 0
            for k in range(N + 1):
                sum += matr[i][k] * self.b[k]
            result.append(sum)
        self.b = result

    def solve(self):
        for i in range(N + 1):
            t = self.make_t(i)
            self.multiply(t)
        return self.b


jordan = JordanGauss(matrix, b)
a = jordan.solve()


def F(x):
    result = 0
    for i in range(N + 1):
        result += a[i] * phi(x, i)
    return result

xF = xFunc
yF = []

for i in range(N + 1):
    yF.append(F(xF[i]))


def Norm(func):
    result = 0
    for i in range(N):
        result += (func(xFunc[i + 1]) ** 2 + func(xFunc[i]) ** 2) * (xFunc[i + 1] - xFunc[i]) / 2.
    return numpy.sqrt(result)


print("Fluff Newton")
print(abs(Norm(f) - Norm(newton_forward)))

print("Fluff Gauss")
print(abs(Norm(f) - Norm(gauss_backward)))

print("Fluff The least squares")
print(abs(Norm(f) - Norm(F)))

# The least square Graphic
plt.plot(xF, yF)

plt.show()