# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import pylab

MIN_X = 0
MAX_X = 1
MIN_Y = -2
MAX_Y = 2
H = 0.125
EXACT_ANSWER_COLOR = "b"
APPROX_ANSWER_COLOR = "r"

N = math.floor((MAX_X - MIN_X) / H)


def exact_answer(x):
    return math.sqrt(x ** 3) - 1


def explicit_euler(y0, grid):
    grid_size = len(grid)
    assert (grid_size >= 2)
    h = grid[1] - grid[0]
    result = [(grid[0], y0)]
    for i in range(grid_size - 1):
        (xi, yi) = result[i]
        yi1 = yi + h * 1.5 * math.sqrt(xi)
        result.append((grid[i + 1], yi1))
    return result


def implicit_euler(y0, grid):
    grid_size = len(grid)
    assert (grid_size >= 2)
    h = grid[1] - grid[0]
    result = [(grid[0], y0)]
    for i in range(grid_size - 1):
        yi = result[i][1]
        xi1 = grid[i + 1]
        yi1 = yi + h * 1.5 * math.sqrt(xi1)
        result.append((xi1, yi1))
    return result


def two_step_adams(y0, grid):
    y1 = (3 * H / 12) * math.sqrt(H) * (1 + 4 * math.sqrt(0.5)) - 1
    result = [(grid[0], y0), (grid[1], y1)]
    for i in range(2, len(grid)):
        y_next = result[i - 1][1] + (H / 2) * (math.sqrt(grid[i - 1]) * 9.0 / 2.0 -
                                               math.sqrt(grid[i - 2]) * 3.0 / 2.0)
        result.append((grid[i], y_next))
    return result


def pairs_to_list(pairs_xy):
    xs, ys = [], []
    for (x, y) in pairs_xy:
        xs.append(x)
        ys.append(y)
    return xs, ys


def draw_axes():
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis([MIN_X, MAX_X, MIN_Y, MAX_Y])
    plt.plot([MIN_X, MAX_X], [0, 0], "k-")


def draw_result(pairs_xy, color):
    xs, ys = pairs_to_list(pairs_xy)
    plt.plot(xs, ys, color + "-")


def generate_grid(x0, x1, n):
    assert (x0 < x1 and n >= 1)
    h = (x1 - x0) / n
    xs = [x0]
    for i in range(n):
        xs.append(x0 + (i + 1) * h)
    return xs


def draw_legend_rect(x, y, color):
    current_axis = plt.gca()
    rect = mpatches.Rectangle(height=0.2, width=0.1, xy=[x, y], alpha=1, color=color)
    current_axis.add_patch(rect)


def test_explicit_euler(x0, x1, n):
    grid = generate_grid(x0, x1, n)
    explicit = explicit_euler(-1.0, grid)
    draw_axes()
    draw_result(explicit, APPROX_ANSWER_COLOR)
    draw_exact_answer(EXACT_ANSWER_COLOR)
    draw_legend_rect(0.54, 1.72, "blue")
    draw_legend_rect(0.54, 1.35, "red")
    pylab.text(0.65, 1.75, u'Точное решение', fontdict={'family': 'verdana'})
    pylab.text(0.65, 1.38, u'Явный метод Эйлера', fontdict={'family': 'verdana'})
    plt.title("h = " + str(H))
    plt.show()


def test_implicit_euler(x0, x1, n):
    grid = generate_grid(x0, x1, n)
    implicit = implicit_euler(-1.0, grid)
    draw_axes()
    draw_result(implicit, APPROX_ANSWER_COLOR)
    draw_exact_answer(EXACT_ANSWER_COLOR)
    draw_legend_rect(0.54, 1.72, "blue")
    draw_legend_rect(0.54, 1.35, "red")
    pylab.text(0.65, 1.75, u'Точное решение', fontdict={'family': 'verdana'})
    pylab.text(0.65, 1.38, u'Неявный метод Эйлера', fontdict={'family': 'verdana'})
    plt.title("h = " + str(H))
    plt.show()


def test_two_step_adams(x0, x1, n):
    grid = generate_grid(x0, x1, n)
    implicit = two_step_adams(-1.0, grid)
    draw_axes()
    draw_result(implicit, APPROX_ANSWER_COLOR)
    draw_exact_answer(EXACT_ANSWER_COLOR)
    draw_legend_rect(0.39, 1.72, "blue")
    draw_legend_rect(0.39, 1.35, "red")
    pylab.text(0.5, 1.75, u'Точное решение', fontdict={'family': 'verdana'})
    pylab.text(0.5, 1.38, u'Явный двушаговый метод Адамса', fontdict={'family': 'verdana'})
    plt.title("h = " + str(H))
    plt.show()


def draw_exact_answer(color):
    xs = generate_grid(MIN_X, MAX_X, 500)
    ys = [exact_answer(x) for x in xs]
    plt.plot(xs, ys, color + "-")

# test_explicit_euler(MIN_X, MAX_X, N)
# test_implicit_euler(MIN_X, MAX_X, N)
test_two_step_adams(MIN_X, MAX_X, N)
