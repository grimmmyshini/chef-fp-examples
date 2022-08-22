#ifndef DUMP_ERRORS_HPP
#define ENABLE_ERROR_DUMP 1

#include <fstream> // Necessary for saving data to plot.
#include <vector>
#include <algorithm>

namespace clad {
    const int no_of_vars = 1;
    const char * var_names[] = {"Ap"};
    const char *output_file_names[] = {"err.Ap.out.dat"};
    std::vector<double> var_errors[no_of_vars];

    void captureVarErrors(const char *name, double error, int iteration = -1) {
        for (int i = 0; i < no_of_vars; i++)
        {
            if (name == var_names[i])
            {
                var_errors[i].push_back(error);
            }
        }
    }

    void dumpVarErrorsToFile() {
        int i = 0, j = 0;
        for (auto &errors : var_errors)
        {
            std::reverse(errors.begin(), errors.end());
            std::ofstream output_file(output_file_names[i]);
            j = 0;
            for (auto &error : errors)
            {
                output_file << j++ << "\t" << error << std::endl;
            }
            i++;
        }
    }
}

#endif