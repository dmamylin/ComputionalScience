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
H = 0.5
EXACT_ANSWER_COLOR = "b"
APPROX_ANSWER_COLOR = "r"

N = math.floor((MAX_X - MIN_X) / H)

def exactAnswer(x):
    return math.sqrt(x**3) - 1

def f(z, xi):
    return z**2 - xi
    
def f_1(z):
    return 2.0 * z

def newton(xi, z0, precision):
    z = z0
    while True:
        zprev = z
        z -= f(z, xi) / f_1(z)
        if abs(z - zprev) < precision:
            break
    return z

def explicitEuler(y0, grid):
    PRECISION = 0.0001
    Z = 1.1
    h = grid[1] - grid[0]
    result = [(grid[0], y0)]
    for i in range(len(grid)):
        yNext = result[i][1] + h * newton(grid[i], Z, PRECISION)
        result.append((grid[i], yNext))
    return result

def implicitEuler(y0, grid):
    PRECISION = 0.0001
    Z = 1.1
    h = grid[1] - grid[0]
    result = [(grid[0], y0)]
    for i in range(1, len(grid)):
        yNext = result[i-1][1] + h * newton(grid[i], Z, PRECISION)
        result.append((grid[i], yNext))
    return result

def twoStepAdams(y0, grid):
    PRECISION = 0.0001
    Z = 1.1
    y1 = (3*H / 12)*math.sqrt(H)*(1 + 4*math.sqrt(0.5)) - 1
    result = [(grid[0], y0), (grid[1], y1)]
    for i in range(2, len(grid)):
        yNext = result[i-1][1] + (H/2)*(newton(grid[i-1], Z, PRECISION)*9.0/2.0 -
                                        newton(grid[i-2], Z, PRECISION)*3.0/2.0)
        result.append((grid[i], yNext))
    return result

def pairsToLists(pairsXY):
    xs, ys = [], []
    for (x, y) in pairsXY:
        xs.append(x)
        ys.append(y)
    return xs, ys

def drawAxes():
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis([MIN_X, MAX_X, MIN_Y, MAX_Y])
    plt.plot([MIN_X, MAX_X], [0, 0], "k-")

def drawResult(pairsXY, color):
    xs, ys = pairsToLists(pairsXY)    
    plt.plot(xs, ys, color+"-")

def generateGrid(x0, x1, n):
    assert(x0 < x1 and n >= 1)
    h = (x1 - x0) / n
    xs = [x0]
    for i in range(n):
        xs.append(x0 + (i+1)*h)
    return xs

def drawLegendRect(x, y, color):
    currentAxis = plt.gca()
    rect = mpatches.Rectangle(height=0.2, width=0.1, xy=[x, y], alpha=1, color=color)
    currentAxis.add_patch(rect)

def testExplicitEuler(x0, x1, n):
    gridX = generateGrid(x0, x1, n)
    explicit = explicitEuler(-1.0, gridX)
    drawAxes()
    drawResult(explicit, APPROX_ANSWER_COLOR)
    drawExactAnswer(EXACT_ANSWER_COLOR)
    drawLegendRect(0.54, 1.72, "blue")
    drawLegendRect(0.54, 1.35, "red")
    pylab.text(0.65, 1.75, u'Точное решение', fontdict={'family':'verdana'})
    pylab.text(0.65, 1.38, u'Явный метод Эйлера', fontdict={'family':'verdana'})
    plt.title("h = " + str(H))
    plt.show()

def testImplicitEuler(x0, x1, n):
    gridX = generateGrid(x0, x1, n)
    implicit = implicitEuler(-1.0, gridX)
    drawAxes()
    drawResult(implicit, APPROX_ANSWER_COLOR)
    drawExactAnswer(EXACT_ANSWER_COLOR)
    drawLegendRect(0.54, 1.72, "blue")
    drawLegendRect(0.54, 1.35, "red")
    pylab.text(0.65, 1.75, u'Точное решение', fontdict={'family':'verdana'})
    pylab.text(0.65, 1.38, u'Неявный метод Эйлера', fontdict={'family':'verdana'})
    plt.title("h = " + str(H))
    plt.show()

def testTwoStepAdams(x0, x1, n):
    gridX = generateGrid(x0, x1, n)
    implicit = twoStepAdams(-1.0, gridX)
    drawAxes()
    drawResult(implicit, APPROX_ANSWER_COLOR)
    drawExactAnswer(EXACT_ANSWER_COLOR)
    drawLegendRect(0.39, 1.72, "blue")
    drawLegendRect(0.39, 1.35, "red")
    pylab.text(0.5, 1.75, u'Точное решение', fontdict={'family':'verdana'})
    pylab.text(0.5, 1.38, u'Явный двушаговый метод Адамса', fontdict={'family':'verdana'})
    plt.title("h = " + str(H))
    plt.show()

def drawExactAnswer(color):
    xs = generateGrid(MIN_X, MAX_X, 500)
    ys = [exactAnswer(x) for x in xs]
    plt.plot(xs, ys, color+"-")    

#testExplicitEuler(MIN_X, MAX_X, N)
#testImplicitEuler(MIN_X, MAX_X, N)
testTwoStepAdams(MIN_X, MAX_X, N)