# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt

PRECISION = 1e-6
MAX_X = 1
MIN_X = 0
H = 0.5
N = math.floor((MAX_X - MIN_X) / H)


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


def f1(x, y, z):
    return z


def f2(x, y, z):
    return y + 16.8 + 7.4 * x * (1.0 - x)


def g1(x, u, v):
    return v


def g2(x, u, v):
    return u


def kutta_solution(grid):
    mu = 1
    grid = generate_grid(MIN_X, MAX_X, N)
    while True:
        shoot = runge_kutta_4(grid, 0, mu, f1, f2)
        return shoot


EXACT_GRAPH = [(x, exact_solution(x)) for x in generate_grid(0, 1, 500)]
KUTTA_GRAPH = [(x, y) for (x, y, z) in kutta_solution(generate_grid(MIN_X, MAX_X, N))]

FUNCTIONS_TO_DRAW = [
    (EXACT_GRAPH, "r"),
    (KUTTA_GRAPH, "r")
]


def draw_axes(xstart, xend, ystart, yend):
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis([xstart, xend, ystart, yend])
    plt.plot([xstart, xend], [0, 0], "k-")
    plt.plot([0, 0], [ystart, yend], "k-")


def draw_function(graph_and_color):
    xs, ys = graph_to_lists(graph_and_color[0])
    plt.plot(xs, ys, graph_and_color[1] + "-")


def draw_functions(list_of_functions):
    draw_axes(0.0, 1.0, -15.0, 15.0)
    for graph_and_color in list_of_functions:
        draw_function(graph_and_color)
    plt.show()


draw_functions(FUNCTIONS_TO_DRAW)
