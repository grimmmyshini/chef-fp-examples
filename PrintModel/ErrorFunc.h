#include <iostream>
#include <limits>

namespace clad {
    double getErrorVal(double dx, double x, std::string name) {
        double error = std::abs(dx * x * std::numeric_limits<float>::epsilon());
        return error;
    }
}
