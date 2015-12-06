from matplotlib import pyplot as plt


def f(x):
    return 2 * x**2


def p(x):
    return -0.5


def q(x):
    return 3.0


a = 1.0
b = 1.3
alpha0 = 1.0
alpha1 = 2.0
beta0 = 1.0
beta1 = 0.0
A = 0.6
B = 1

def solve(n):
    xval = []
    yval = []
    c = []
    d = []
    h = (b - a) / n
    c.append(alpha1 / (h * alpha0 - alpha1))
    d.append(A * h / (h * alpha0 - alpha1))
    xval.append(a)
    for i in range(1, n):
        xval.append(a + i * h)
        mi = ((h ** 2) * q(xval[i]) - 2) / (1 + h * p(xval[i]) / 2)
        ki = (1 - h * p(xval[i]) / 2) / (1 + h * p(xval[i]) / 2)
        Fi = f(xval[i]) / (1 + h * p(xval[i]) / 2)
        c.append(1 / (mi - ki * c[i - 1]))
        d.append(((h ** 2) * Fi - ki * d[i - 1]) * c[i])
    xval.append(b)
    yval.append((B * h + beta1 * d[n - 1]) / (beta0 * h + beta1 * (c[n - 1] + 1)))
    for i in range(n):
        yval.append(d[n - 1 - i] - c[n - 1 - i] * yval[i])
    yval.reverse()
    return xval, yval

xcoords, ycoords = solve(100)
plt.plot(xcoords, ycoords)
plt.show()