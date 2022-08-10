#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>

namespace clad
{

    class ErrorStorage
    {
    private:
        struct Properties
        {
            double total_error;
            long long count;
            double max_error;
        };

    public:
        static ErrorStorage &getInstance()
        {
            static ErrorStorage instance;
            return instance;
        }

        __attribute__((always_inline)) void set_error(const char *var, const double err)
        {
            auto entry = error_map.find(var);

            if (entry != error_map.end())
            {
                entry->second.total_error += err;
                assert(entry->second.count > 0 && "Count exceeds limits!");
                entry->second.count++;
                if (entry->second.max_error < err)
                    entry->second.max_error = err;
            }
            else
            {
                error_map[var] = {err, 1, err};
            }
        }

        __attribute__((always_inline)) void reset()
        {
            error_map.clear();
        }

        const std::unordered_map<const char *, Properties> &get_error_map() const
        {
            return error_map;
        }

        ErrorStorage(ErrorStorage const &) = delete;
        void operator=(ErrorStorage const &) = delete;

    private:
        ErrorStorage() {}

        std::unordered_map<const char *, Properties> error_map;

        friend class ErrorReport;
    };

    class ErrorReport
    {
    public:
        ErrorReport(std::unordered_map<const char *, ErrorStorage::Properties> error_map)
            : reports(), max_var_name_size(0)
        {
            for (const auto &error : error_map)
            {
                std::string var_name(error.first);
                double total_error = error.second.total_error;
                long long count = error.second.count;
                double max_error = error.second.max_error;
                double avg_error;
                bool to_replace;
                max_var_name_size = max_var_name_size < var_name.size() ? var_name.size() : max_var_name_size;

                if (count > 0)
                    avg_error = total_error / count;
                else
                    avg_error = std::numeric_limits<double>::quiet_NaN();

                to_replace = (total_error < threshold);

                reports.push_back({var_name, total_error, count, max_error, avg_error, to_replace});
            }
        }

        void print()
        {
            std::sort(reports.begin(), reports.end(), compareByTotalError);
            std::sort(reports.begin(), reports.end(), compareByToReplace);

            max_var_name_size = max_var_name_size > 8 ? max_var_name_size : 8;

            std::string margin("  ");

            long long length = 7 + max_var_name_size + 11 + 11 + 11 + 11 + margin.size() * 5;

            print_banner(" CLAD ERROR REPORT ", length, '=');

            print_with_padding("Replace", 7);
            std::cout << margin;
            print_with_padding("Var Name", max_var_name_size);
            std::cout << margin;
            print_with_padding("Max Error", 11);
            std::cout << margin;
            print_with_padding("Count", 11);
            std::cout << margin;
            print_with_padding("Total Error", 11);
            std::cout << margin;
            print_with_padding("Avg. Error", 11);
            std::cout << std::endl;
            print_banner("", length, '-');

            for (auto &report : reports)
            {
                print_with_padding(report.to_replace ? "YES" : "NO", 7);
                std::cout << margin;
                print_with_padding(report.var_name, max_var_name_size);
                std::cout << margin;
                print_with_padding(report.max_error, 11);
                std::cout << margin;
                print_with_padding(report.count, 11);
                std::cout << margin;
                print_with_padding(report.total_error, 11);
                std::cout << margin;
                print_with_padding(report.avg_error, 11);
                std::cout << std::endl;
            }

            print_banner("", length, '=');
        }

    private:
        constexpr const static double threshold = 1e-12;

        struct Report
        {
            std::string var_name;
            double total_error;
            long long count;
            double max_error;
            double avg_error;
            bool to_replace;
        };

        std::vector<Report> reports;
        long long max_var_name_size;

        static bool compareByTotalError(const Report &a, const Report &b)
        {
            return a.total_error < b.total_error;
        }

        static bool compareByToReplace(const Report &a, const Report &b)
        {
            return a.to_replace > b.to_replace;
        }

        static void print_with_padding(std::string str, long long length)
        {
            std::cout << str;

            for (long long i = 0, end = length - str.size(); i < end; i++)
                std::cout << " ";
        }

        static void print_with_padding(double num, long long length)
        {
            std::cout.precision(length - 6);
            std::cout << std::scientific << num;
        }

        static void print_with_padding(long long num, long long length)
        {
            auto num_str = std::to_string(num);
            print_with_padding(num_str, length);
        }

        static void print_banner(std::string str, long long length, char symbol)
        {
            for (long long i = 0, end = (length - str.size()) / 2; i < end; i++)
                std::cout << symbol;
            std::cout << str;
            for (long long i = 0, end = (length - str.size() + 1) / 2; i < end; i++)
                std::cout << symbol;
            std::cout << std::endl;
        }
    };

    void printErrorReport()
    {
        auto error_map = ErrorStorage::getInstance().get_error_map();
        ErrorReport error_report(error_map);

        error_report.print();
    }

    void resetErrors()
    {
        ErrorStorage::getInstance().reset();
    }

    __attribute__((always_inline)) double getErrorVal(double dx, double x, const char *name)
    {
        double error = std::abs(dx * x * std::numeric_limits<float>::epsilon());
        ErrorStorage::getInstance().set_error(name, error);
        return error;
    }
} // clad
