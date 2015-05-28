#include <stdio.h>
#include <math.h>

#define A 2.0
#define B 3.0

double f(double x) {
    return exp(-cos(x));
}

double midRectangles(double a, double b, double h) {
    double result = 0.0;
    int i;
    int n = (int)floor((b - a) / h);

    for (i = 0; i < n; i++) {
        double xPrev = a + i * h;
        double xNext = a + (i + 1) * h;
        double middlePoint = (xPrev + xNext) / 2.0;

        result += f(middlePoint);
    }

    return result * h;
}

double gregory(double a, double b, double h) {
    double x1 = a + h;
    double xn_1 = b - h;
    double result = 0.0;
    int i;
    int n = (int)floor((b - a) / h);

    result += (5 * (f(a) + f(b)) + f(x1) + f(xn_1)) / 12.0;
    
    for (i = 1; i < n; i++) {
        result += f(a + i * h);
    }

    return result * h;
}

int main(int argc, char** argv) {
    int i;
    const int stepsCount = 3;
    const double steps[] = { 0.1, 0.05, 0.025 };

    for (i = 0; i < stepsCount; i++) {
        printf("Step = %.3lf:\n", steps[i]);
        printf("\tRigth triangles = %.6lf\n", midRectangles(A, B, steps[i]));
        printf("\tGregory = %.6lf\n", gregory(A, B, steps[i]));
    }

    return 0;
}