#include <cmath>
#include "adapt.h"

#define ABS(x) (((x) < 0.0) ? (-(x)) : (x))

#define PI 3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

#define ITERATIONS 1000000

namespace adapt
{

    AD_real h = PI / (double)ITERATIONS; // double

    AD_real t1 = 0.0;
    AD_real t2;
    AD_real t3; // double

    AD_real s1 = 0.0; // double

    AD_real d1 = 1.0;
    AD_real d2; // float

    AD_real fun(AD_real x)
    {
        d2 = d1; // also d1 in original
        AD_INTERMEDIATE(d2, "d2");
        t3 = x; // also t1 in original
        AD_INTERMEDIATE(t3, "t3");

        int k;
        for (k = 1; k <= 5; k += 1)
        {
            d2 = 2.0 * d2;
            AD_INTERMEDIATE(d2, "d2");
            t3 = t3 + sin(d2 * x) / d2;
            AD_INTERMEDIATE(t3, "t3");
        }
        return t3;
    }

    void do_fun()
    {
        int i;
        for (i = 1; i <= ITERATIONS; i += 1)
        {
            t2 = fun(i * h);
            AD_INTERMEDIATE(t2, "t2");
            s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
            AD_INTERMEDIATE(s1, "s1");
            t1 = t2;
            AD_INTERMEDIATE(t1, "t1");
        }
    }

} // adapt

namespace clad
{

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

}