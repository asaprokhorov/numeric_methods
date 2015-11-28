from matplotlib import pyplot as plt
import math

a = 0.29
b = 1.29
y0 = 0.34
eps = 0.00001


def f(x, y):
    return 0.158 * (x ** 2 + math.cos(0.6 * x)) + 0.675 * y


def eiler(n):
    xval = []
    yval = []
    h = (b - a) / n
    xval.append(a)
    yval.append(y0)
    for i in range(1, n + 1):
        xval.append(a + i * h)
        yi = yval[i - 1] + h * f(xval[i - 1], yval[i - 1])
        yi1 = yval[i - 1] + h * (f(xval[i - 1], yval[i - 1]) + f(xval[i - 1], yi)) / 2
        while(abs(yi1 - yi) >= eps):
            yi = yi1
            yi1 = yi + h * (f(xval[i - 1], yval[i - 1]) + f(xval[i - 1], yi)) / 2
        yval.append(yi1)
    return xval, yval


def rungeKutta(n):
    xval = []
    yval = []
    h = (b - a) / n
    xval.append(a)
    yval.append(y0)
    for i in range(1, n + 1):
        xval.append(a + i * h)
        k1 = h * f(xval[i - 1], yval[i - 1])
        k2 = h * f(xval[i - 1] + (1/2) * h, yval[i - 1] + (1/2) * k1)
        k3 = h * f(xval[i - 1] + (1/2) * h, yval[i - 1] + (1/2) * k2)
        k4 = h * f(xval[i - 1] + (1/2) * h, yval[i - 1] + k3)
        yn = yval[i - 1] + (1/6) * (k1 + 2 *  k2 + 2 * k3 + k4)
        yval.append(yn)
    return xval, yval


def gauss(s, funct):
    result = 0
    t = [0.906, 0.538469, 0, -0.538469, -0.90618]
    C = [0.23693, 0.47863, 0.56889, 0.47863, 0.23693]
    for i in range(0, 5):
        result += C[i] * funct(s, (1.0 + 0.0) / 2 + (1.0 - 0.0) / 2 * t[i])
    result *= (1.0 - 0.0) / 2
    return result


def differences(n):
    xval, yval = rungeKutta(n)
    result = []
    funct = []
    for i in range(n + 1):
        funct.append(f(xval[i], yval[i]))
    result.append(funct);
    for i in range(1, n + 1):
        result.append([])
        for j in range(0, n - i + 1):
            result[i].append(result[i - 1][j + 1] - result[i - 1][j])
    return result


def alphaE(s, t):
    result = 1.0
    if s == 0:
        return result
    for i in range(0, s):
       result *= t + i
    result /= math.factorial(s)
    return result


def extrapolationAdams(n, m):
    diffs = differences(n)
    xval, yval = rungeKutta(n)
    h = (b - a) / n
    for i in range(m + 1, n + 1):
        sum = 0
        for s in range(m + 1):
            sum += gauss(s, alphaE) * diffs[s][i - 1 - s]
        yval[i] = yval[i - 1] + h * sum
    return xval, yval


def alphaI(s, t):
    result = 1.0
    for i in range(0, s):
        result *= t + i - 1
    result /= math.factorial(s)
    return result


def interpolationAdams(n, m):
    diffs = differences(n)
    xval, yval = rungeKutta(n)
    h = (b - a) / n
    for i in range(m + 1, n + 1):
        sum = 0
        for s in range(1, m + 1):
            sum += gauss(s, alphaI) * diffs[s][i - s]
        d = yval[i - 1] + h * sum
        yi = yval[i - 1]
        yi1 = h * f(xval[i], yi) + d
        while abs(yi1 - yi) >= eps:
            yi = yi1
            yi1 = h * f(xval[i], yi) + d
        yval[i] = yi1
    return xval, yval


x_eiler, y_eiler = eiler(100)
x_runge, y_runge = rungeKutta(100)
x_extrapolation, y_extrapolation = extrapolationAdams(100, 49)
x_interpolation, y_interpolation = interpolationAdams(100, 49)

plt.plot(x_eiler, y_eiler)
plt.plot(x_runge, y_runge)
plt.plot(x_extrapolation, y_extrapolation)
plt.plot(x_interpolation, y_interpolation)
plt.show()