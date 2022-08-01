#include <iostream>

namespace clad {
    double getErrorVal(double dx, double x, std::string name) {
        double error = std::abs(dx * (x - (float)x));
        std::cout << name << " : " << error << std::endl;
        return error;
    }
}