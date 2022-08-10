#ifndef ARCLEN_HPP
#define ARCLEN_HPP

#include <cmath>
#define ABS(x) (((x) < 0.0) ? (-(x)) : (x))

#define PI 3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

#define ITERATIONS 1000000

double do_fun(double s1, double t1)
{
    int i, k;
    double t2, h = PI / ITERATIONS, x, d2, t3;

    for (i = 1; i <= ITERATIONS; i += 1)
    {
        // t2 = fun_clad(i * h);
        x = i * h;
        d2 = 1.0; // also d1 in original
        t3 = x;   // also t1 in original

        for (k = 1; k <= 5; k += 1)
        {
            d2 = 2.0 * d2;
            t3 = t3 + sin(d2 * x) / d2;
        }

        t2 = t3;

        s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
        t1 = t2;
    }

    return s1;
}

#endif