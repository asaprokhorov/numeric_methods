import math
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


print("ITERATIONSYS")
print(iterationsys(0, 0, 0.00001))

print("\nNEWTONSYS")
print(newtonsys(0, 0, 0.00001))

