#include <iomanip>  // necessary for setprecision
#include <iostream> // necessary for 'cout'

#include "simpsons-adapt.hpp"

#include "adapt.h"

int main() {
  AD_begin();

  AD_real a = 0, b = 1;

  AD_INDEPENDENT(a, "a");
  AD_INDEPENDENT(b, "b");
  
  AD_real result = adapt::simpsons(a, b, 1000);

  AD_DEPENDENT(result, "result", 0.0);
  AD_report();
  AD_end();
}
