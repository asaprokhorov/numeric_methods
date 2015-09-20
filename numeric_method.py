import math


def f(x):
    return x * x + 4 * math.sin(x)


def prime(x):
    return 2 * x + 4 * math.cos(x)


def chords(a, b):
    return (a * f(b) - b * f(a)) / (f(b) - f(a))


def newton(a, b):
    return b - f(b)/prime(b)


def combined(a, b):
    x = newton(a, b)
    return chords(x, b)


def linear(a, b, eps, method):
    x = method(a, b)
    iteration = 0
    while abs(x - b) > eps:
        b = x
        iteration += 1
        x = method(a, b)
    return x, iteration


def s1(x, y):
    return math.sin(x - 1) - 1.3 + y


def s2(x, y):
    return x - math.sin(y + 1) - 0.8


def s1x(x, y):
    return math.cos(x - 1)


def s1y(x, y):
    return 1


def s2x(x, y):
    return 1


def s2y(x, y):
    return math.cos(y + 1)


def dn(x, y):
    return s1x(x, y) * s2y(x, y) - s1y(x, y) * s2x(x, y)


def dx(x, y):
    return s1y(x, y) * s2(x, y) - s1(x, y) * s2y(x, y)


def dy(x, y):
    return s1(x, y) * s2x(x, y) - s1x(x, y) * s2(x, y)


def newtonsys(x0, y0, eps):
    x = x0 + dx(x0, y0) / dn(x0, y0)
    y = y0 + dy(x0, y0) / dn(x0, y0)
    iteration = 0
    while abs(x - x0) > eps and abs(y - y0) > eps:
        x0 = x
        y0 = y
        iteration += 1
        x = x0 + dx(x0, y0) / dn(x0, y0)
        y = y0 + dy(x0, y0) / dn(x0, y0)
    return x, y, iteration


def itY(x, y):
    return 1.3 - math.sin(x - 1)


def itX(x, y):
    return math.sin(y + 1) + 0.8


def iterationsys(x0, y0, eps):
    x = itX(x0, y0)
    y = itY(x0, y0)
    iteration = 0
    while abs(x - x0) > eps and abs(y - y0) > eps:
        iteration += 1
        x0 = x
        y0 = y
        x = itX(x0, y0)
        y = itY(x0, y0)
    return x, y, iteration


print("CHORDS")
print(linear(-0.05, 0.05, 0.00000001, chords))
print("NEWTON")
print(linear(-0.05, 0.05, 0.00000001, newton))
print("COMBINED")
print(linear(-0.05, 0.05, 0.00000001, combined))

print("ITERATIONSYS")
print(iterationsys(0, 0, 0.00001))

print("\nNEWTONSYS")
print(newtonsys(0, 0, 0.00001))

