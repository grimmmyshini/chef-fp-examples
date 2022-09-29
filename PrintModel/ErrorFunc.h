#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <utility>
#include <cassert>
#include <algorithm>

// Fast approx
#include "fastonebigheader.h"

namespace clad
{
    void captureVarErrors(const char *name, double error, int iteration = -1);
    void dumpVarErrorsToFile();
    
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
                bool to_lower_precision;
                max_var_name_size = max_var_name_size < var_name.size() ? var_name.size() : max_var_name_size;

                if (count > 0)
                    avg_error = total_error / count;
                else
                    avg_error = std::numeric_limits<double>::quiet_NaN();

                to_lower_precision = (total_error < threshold);

                reports.push_back({var_name, total_error, count,
                                   max_error, avg_error, to_lower_precision});
            }
        }

        void print()
        {
            std::sort(reports.begin(), reports.end(), compareByTotalError);
            std::sort(reports.begin(), reports.end(), compareByToLowerPrecision);

            max_var_name_size = max_var_name_size > 8 ? max_var_name_size : 8;

            std::string margin("  ");

            const int len_prec = 9;
            const int len_varn = max_var_name_size;
            const int len_merr = 11;
            const int len_coun = 11;
            const int len_terr = 11;
            const int len_aerr = 11;

            long long length = len_prec + len_varn + len_merr + len_coun +
                               len_terr + len_aerr + margin.size() * 5;

            print_banner(" CLAD ERROR REPORT ", length, '=');

            print_with_padding("Precision", len_prec);
            std::cout << margin;
            print_with_padding("Var Name", len_varn);
            std::cout << margin;
            print_with_padding("Max Error", len_merr);
            std::cout << margin;
            print_with_padding("Count", len_coun);
            std::cout << margin;
            print_with_padding("Total Error", len_terr);
            std::cout << margin;
            print_with_padding("Avg. Error", len_aerr);
            std::cout << std::endl;
            print_banner("", length, '-');

            for (auto &report : reports)
            {
                print_with_padding(report.to_lower_precision ? "lower" : "HIGHER", len_prec);
                std::cout << margin;
                print_with_padding(report.var_name, len_varn);
                std::cout << margin;
                print_with_padding(report.max_error, len_merr);
                std::cout << margin;
                print_with_padding(report.count, len_coun);
                std::cout << margin;
                print_with_padding(report.total_error, len_terr);
                std::cout << margin;
                print_with_padding(report.avg_error, len_aerr);
                std::cout << std::endl;
            }

            print_banner("", length, '=');
        }

    private:
        constexpr const static double threshold = 1e-6;

        struct Report
        {
            std::string var_name;
            double total_error;
            long long count;
            double max_error;
            double avg_error;
            bool to_lower_precision;
        };

        std::vector<Report> reports;
        long long max_var_name_size;

        static bool compareByTotalError(const Report &a, const Report &b)
        {
            return a.total_error < b.total_error;
        }

        static bool compareByToLowerPrecision(const Report &a, const Report &b)
        {
            return a.to_lower_precision > b.to_lower_precision;
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

        #ifdef ENABLE_ERROR_DUMP
        dumpVarErrorsToFile();
        #endif
    }

    void resetErrors()
    {
        ErrorStorage::getInstance().reset();
    }

    __attribute__((always_inline)) double getErrorVal(double dx, double x, const char *name)
    {
        double error = std::abs(dx * (x - (float)x)); // replace x * std::numeric_limits<float>::epsilon() with avg or max error
        ErrorStorage::getInstance().set_error(name, error);
        #ifdef ENABLE_ERROR_DUMP
        captureVarErrors(name, error);
        #endif
        return error;
    }

    __attribute__((always_inline)) double doApprox(double dx, double x, const char* cname) {
      char name[50];
      int i = 0;
      while (cname[i++] != '\0')
        name[i - 1] = cname[i - 1];
      char *token = strtok(name, "_");
      if (strcmp(token, "clad"))
        return 0;
      token = strtok(NULL, "_");
      double error;
      if (!strcmp(token, "exp"))
        error = std::fabs(dx * (exp(x) - fastexp(x)));
      else if (!strcmp(token, "log"))
        error = std::fabs(dx * (log(x) - fastlog(x)));
      else if (!strcmp(token, "sqr"))
        error = std::fabs(dx * (sqrt(x) - fastpow(x, 0.5)));
      else
        return 0;

      ErrorStorage::getInstance().set_error(cname, error);
      return error;
    }
} // clad
