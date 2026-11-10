/*
 * Вариант 16. Задача Дирихле для уравнения Пуассона
 *     -Δu = f   в D,      u = φ   на γ = ∂D,
 * D = {x^2 + y^2 < 2} ∩ {|y| < 1}  (часть круга внутри полосы),
 * f(x,y) = 1,  φ(x,y) = sin(x + y)
 *
 * Метод фиктивных областей: D вкладывается в прямоугольник
 *     Π = [-√2, √2] × [-1, 1],
 * в фиктивной области Π \ D ставится Q = 1/ε^2, F = 0, ε = h^2
 *
 * При компиляции с -DCHECK_EXACT решается тестовая задача с известным
 * решением u = x^2 + y^2 (f = -4, φ = x^2 + y^2) — удобно для проверки кода
 */
#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>

#define SQRT2 1.41421356237309504880

// Границы объемлющего прямоугольника Π = [A1, B1] × [A2, B2]
static const double A1 = -SQRT2, B1 = SQRT2;
static const double A2 = -1.0,   B2 = 1.0;

#ifdef CHECK_EXACT
static double f_rhs(double x, double y)
{
    (void)x;
    (void)y;
    return -4.0;
}

static double phi(double x, double y)
{
    return x * x + y * y;
}
#else
static double f_rhs(double x, double y)
{
    (void)x;
    (void)y;
    return 1.0;
}

static double phi(double x, double y)
{
    return sin(x + y);
}
#endif

// Принадлежит ли точка (x, y) исходной области D
static int in_domain(double x, double y)
{
    return x * x + y * y < 2.0 && fabs(y) < 1.0;
}

// Длина пересечения вертикального отрезка {x} × [ya, yb] с областью D
static double seg_len(double x, double ya, double yb)
{
    double g, lo, hi;

    // Вне круга по x — сечение пустое
    if (x * x >= 2.0) {
        return 0.0;
    }

    // Верхняя граница по y: минимум из полосы |y| < 1 и окружности
    g = sqrt(2.0 - x * x);
    if (g > 1.0) {
        g = 1.0;
    }

    // Пересекаем [ya, yb] с симметричным интервалом [-g, g]
    lo = ya > -g ? ya : -g;
    hi = yb <  g ? yb :  g;
    return hi > lo ? hi - lo : 0.0;
}

/*
 * Площадь пересечения ячейки [xa, xb] × [ya, yb] с областью D
 * Область выпукла, поэтому если все четыре угла лежат в D, ячейка
 * целиком внутри и её площадь известна точно. Иначе интегрируем длину
 * вертикального сечения по x составной формулой средних (по y — точно)
 */
static double cell_area(double xa, double xb, double ya, double yb)
{
    const int K = 128;
    double dx, s = 0.0;
    int k;

    if (in_domain(xa, ya) && in_domain(xa, yb) &&
        in_domain(xb, ya) && in_domain(xb, yb)) {
        return (xb - xa) * (yb - ya);
    }

    dx = (xb - xa) / K;
    for (k = 0; k < K; k++) {
        s += seg_len(xa + (k + 0.5) * dx, ya, yb);
    }
    return s * dx;
}

#endif
