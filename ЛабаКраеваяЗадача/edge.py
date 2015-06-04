# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab

PRECISION = 1e-6
MAX_X = 1
MIN_X = 0
H = 0.05
N = math.floor((MAX_X - MIN_X) / H)


def f1(x, y, z):
    return z


def f2(x, y, z):
    return y + 16.8 + 7.4 * x * (1.0 - x)


def g1(x, u, v):
    return v


def g2(x, u, v):
    return u

def exact_solution(x):
    double_e = 2.0 * math.e
    return math.exp(-x)*(double_e + 3.7) / double_e +\
           math.exp(x)*(double_e - 3.7) / double_e +\
           7.4*(x**2 - x) - 2


def generate_grid(a, b, n):
    h = (b - a) / n
    return [a + i * h for i in range(n + 1)]


def graph_to_lists(graph):
    xs, ys = [], []
    for (x, y) in graph:
        xs.append(x)
        ys.append(y)
    return xs, ys


"""Метод Рунге-Кутта четвертого порядка для системы из двух ОДУ"""
def runge_kutta_4(grid, y0, z0, f1, f2):
    n = len(grid)
    assert (n >= 2)
    h = grid[1] - grid[0]
    result = [(grid[0], y0, z0)]
    for i in range(1, n):
        xprev = result[i - 1][0]
        yprev = result[i - 1][1]
        zprev = result[i - 1][2]
        k1 = h * f1(xprev, yprev, zprev)
        l1 = h * f2(xprev, yprev, zprev)

        k2 = h * f1(xprev + h / 2.0, yprev + k1 / 2.0, zprev + l1 / 2.0)
        l2 = h * f2(xprev + h / 2.0, yprev + k1 / 2.0, zprev + l1 / 2.0)

        k3 = h * f1(xprev + h / 2.0, yprev + k2 / 2.0, zprev + l2 / 2.0)
        l3 = h * f2(xprev + h / 2.0, yprev + k2 / 2.0, zprev + l2 / 2.0)

        k4 = h * f1(grid[i], yprev + k3, zprev + l3)
        l4 = h * f2(grid[i], yprev + k3, zprev + l3)

        ynext = yprev + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        znext = zprev + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
        result.append((grid[i], ynext, znext))
    return result


def newton(x0, f, f_derivative, precision):
    result = x0
    while True:
        xprev = result
        result = xprev - f(xprev) / f_derivative(xprev)
        if math.fabs(result - xprev) < precision:
            break
    return result


"""Применение метода стрельбы"""
def shoot_method(grid):
    mu = 1
    while True:
        shoot = runge_kutta_4(grid, 0, mu, f1, f2)
        mu_shoot = runge_kutta_4(grid, 0, 1, g1, g2)
        n = len(shoot)

        # Находим F и F':
        big_f = shoot[n-1][1] + shoot[n-1][2] - 2.0*math.e - 1.7
        big_f_derivative = mu_shoot[n-1][1] + mu_shoot[n-1][2]

        # Метод Ньютона:
        new_mu = mu - big_f/big_f_derivative

        if math.fabs(new_mu - mu) < PRECISION:
            return shoot
        mu = new_mu


ABOVE_MAIN_DIAGONAL = [1.0 for x in range(N)]
ABOVE_MAIN_DIAGONAL[0] = 0.0

BELOW_MAIN_DIAGONAL = [1.0 for x in range(N)]
BELOW_MAIN_DIAGONAL[N-1] = 2.0

MAIN_DIAGONAL = [-H**2 - 2.0 for x in range(N+1)]
MAIN_DIAGONAL[N] = -2.0 - H**2 - 2.0*H
MAIN_DIAGONAL[0] = 1.0

RIGHT_PART = [(16.8 + 7.4*xi*(1 - xi))*H**2 for xi in generate_grid(MIN_X, MAX_X, N)]
RIGHT_PART[0] = 0.0
RIGHT_PART[N] = 16.8*H**2 - 2.0*H*(2.0*math.e + 1.7)
"""def tridiagonal_matrix_method(below, above, main, right, grid):
    y0 = (right[0] - above[0]) / main[0]
    y1 = right[1] / ABOVE_MAIN_DIAGONAL[0]
    n = len(main)
    result = [y0, y1]

    for i in range(1, n-1):
        yi = (right[i] - below[i-1]*result[i-1] - main[i]*result[i]) / above[i]
        result.append(yi)

    # yn = (16.8*H**2 - 2*math.e - 1.7 - 2.0*result[n-2])/(-2.0 - H**2 - 2.0*H)
    yn = (right[n-1] - below[n-2]*result[n-1]) / main[n-1]
    result.append(yn)

    return zip(grid, result)"""
def tridiagonal_matrix_method(below, above, main, right, grid):
    a, b, c, d = below, main, above, right
    n = len(main)

    a[0] /= b[0]
    for i in range(1, n-1):
        b[i] -= (a[i-1]*c[i-1])
        a[i] /= b[i]

    b[n-1] -= a[n-2]*c[n-2]

    for i in range(1, n):
        d[i] -= d[i-1]*a[i-1]
    d[n-1] /= b[n-1]

    for i in range(n-2, -1, -1):
        d[i] -= d[i+1]*c[i]
        d[i] /= b[i]

    return [(grid[i], d[i]) for i in range(len(grid))]


EXACT_GRAPH = [(x, exact_solution(x)) for x in generate_grid(0, 1, 500)]
SHOOT_GRAPH = [(x, y) for (x, y, z) in shoot_method(generate_grid(MIN_X, MAX_X, N))]
DIAG_GRAPH = tridiagonal_matrix_method(BELOW_MAIN_DIAGONAL, ABOVE_MAIN_DIAGONAL, MAIN_DIAGONAL, RIGHT_PART,
                                       generate_grid(MIN_X, MAX_X, N))

FUNCTIONS_TO_DRAW = [
    (EXACT_GRAPH, "r"),
    (DIAG_GRAPH, "b"),
    (SHOOT_GRAPH, "g")
]

LEGENDS_TO_DRAW = [
    (0.2, -0.2, "r", u"Точное решение"),
    (0.2, -0.4, "g", u"Метод стрельбы"),
    (0.2, -0.6, "b", u"Метод прогонки")
]


def draw_legend_rect(x, y, color):
    current_axis = plt.gca()
    rect = mpatches.Rectangle(height=0.15, width=0.05, xy=[x, y], alpha=1, color=color)
    current_axis.add_patch(rect)


def draw_legend(x, y, color, text):
    draw_legend_rect(x, y, color)
    pylab.text(x + 0.061, y + 0.05, text, fontdict={'family': 'verdana'})


def draw_axes(xstart, xend, ystart, yend):
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis([xstart, xend, ystart, yend])
    plt.plot([xstart, xend], [0, 0], "k-")
    plt.plot([0, 0], [ystart, yend], "k-")


def draw_function(graph_and_color):
    xs, ys = graph_to_lists(graph_and_color[0])
    plt.plot(xs, ys, graph_and_color[1] + "-")


def draw_functions(list_of_functions, list_of_legends, title):
    draw_axes(0.0, 1.0, -2.5, 0.0)
    for graph_and_color in list_of_functions:
        draw_function(graph_and_color)
    for (x, y, color, text) in list_of_legends:
        draw_legend(x, y, color, text)
    plt.title(title)
    plt.show()


draw_functions(FUNCTIONS_TO_DRAW, LEGENDS_TO_DRAW, "h = " + str(H))
