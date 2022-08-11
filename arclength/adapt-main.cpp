#define ITERATIONS 1000000

#include "arclen.hpp"
#include "arclen-adapt.hpp"

#include "adapt.h"



int main()
{
    using namespace adapt;
    AD_begin();
    AD_INTERMEDIATE(h, "h");
    AD_INTERMEDIATE(t1, "t1");
    AD_INTERMEDIATE(s1, "s1");
    AD_INTERMEDIATE(d1, "d1");

    do_fun();

    AD_DEPENDENT(s1, "s1", EPS);
    AD_report();
    t1 = 0.0;
    s1 = 0.0;
    d1 = 1.0;
}