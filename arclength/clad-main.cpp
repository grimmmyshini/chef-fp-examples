#include "arclen.hpp"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "Derivative.hpp"

// #include <gperftools/profiler.h>

int main()
{
    //   auto df = clad::estimate_error(do_fun);
    double ds1 = 0, dt1 = 0;
    double fin_error = 0;

    // ProfilerStart("arclength.prof");
    double result = do_fun(0, 0);
    // ProfilerStop();


    // ProfilerStart("clad.prof");
    clad::resetErrors();
    clad::do_fun_grad(0, 0, &ds1, &dt1, fin_error);
    clad::printErrorReport();
    // ProfilerStop();
}