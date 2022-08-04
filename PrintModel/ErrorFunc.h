#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>

namespace clad
{
    class MaxError
    {
    public:
        static MaxError &getInstance()
        {
            static MaxError instance;
            return instance;
        }

        __attribute__((always_inline)) void set_error(const char * var, const double err)
        {
            auto entry = error_map.find(var);

            if (entry != error_map.end())
            {
                entry->second.first += err;
                assert(entry->second.second > 0 && "Count exceeds limits!");
                entry->second.second++;
            }
            else
            {
                error_map[var] = {err, 1};
            }
        }

        __attribute__((always_inline)) void reset()
        {
            error_map.clear();
        }

        const std::unordered_map<const char *, std::pair<double, long long>> &get_error_map() const
        {
            return error_map;
        }

        MaxError(MaxError const &) = delete;
        void operator=(MaxError const &) = delete;

    private:
        MaxError() {}
        std::unordered_map<const char *, std::pair<double, long long>> error_map;
    };

    void printErrorReport()
    {
        auto error_map = MaxError::getInstance().get_error_map();

        std::cout << "Clad error report: " << std::endl;
        for (const auto &error : error_map)
        {
            double total = error.second.first;
            long long count = error.second.second;

            if (count > 0)
                std::cout << error.first << ":\ttotal = " << total << "\tcount = " << count << "\tavg = " << total / count << std::endl;
            else

                std::cout << error.first << ":\ttotal = " << total << "\tcount = out of bounds"
                          << "\tavg = undeterminable" << std::endl;
        }
    }

    void resetErrors()
    {
        MaxError::getInstance().reset();
    }

    __attribute__((always_inline)) double getErrorVal(double dx, double x, const char* name)
    {
        double error = std::abs(dx * x * std::numeric_limits<float>::epsilon());
        MaxError::getInstance().set_error(name, error);
        return error;
    }
} // clad
