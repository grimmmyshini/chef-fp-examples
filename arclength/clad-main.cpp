#include "arclen.hpp"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "Derivative.hpp"

// #include <gperftools/profiler.h>

int main()
{
    // auto df = clad::estimate_error(do_fun<double>);
    int diter = 0;
    double fin_error = 0;

    // ProfilerStart("arclength.prof");
    // long double ans = 5.795776322412856;
    // std::cout << "Actual Errors: "
    //           << "\nFloat: "
    //           << std::fabs(ans - do_fun<float>(1000000))
    //           << "\nDouble: "
    //           << std::fabs(ans - do_fun<double>(1000000))
    //           << "\nMixed precision: "
    //           << std::fabs(ans - do_fun<double, float>(1000000))
    //           << std::endl;
    // ProfilerStop();

    // ProfilerStart("clad.prof");
    clad::resetErrors();
    clad::do_fun_grad(1000000, &diter, fin_error);
    clad::printErrorReport();
    // std::cout << "Final Estimated Error: " << fin_error << std::endl;
    // ProfilerStop();
}
