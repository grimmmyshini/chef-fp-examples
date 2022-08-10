#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "benchmark/benchmark.h"

#include "generate_matrix.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"

#undef DEBUG

#define nx 30
#define ny 40
#define nz 1600

struct cout_suppressor
{
    cout_suppressor()
        : buffer(), old(std::cout.rdbuf(buffer.rdbuf()))
    {
    }

    ~cout_suppressor()
    {
        std::cout.rdbuf(old);
    }

private:
    std::stringstream buffer;
    std::streambuf *old;
};

HPC_Sparse_Matrix A;

namespace lower_prec
{
    float HPCCG_residual(double *b, double *x, int max_iter, float tolerance,
                         int &niters, float &normr, double *r, float *p,
                         float *Ap, double *xexact)
    {
        // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
        int nrow = A.local_nrow;
        int ncol = A.local_ncol;

        normr = 0.0;
        double rtrans = 0.0;
        float oldrtrans = 0.0;
        double alp, bta, alpha, beta;

        int rank = 0; // Serial case (not using MPI)

        int print_freq = max_iter / 10;
        if (print_freq > 50)
            print_freq = 50;
        if (print_freq < 1)
            print_freq = 1;

        // // p is of length ncols, copy x to p for sparse MV operation
        // waxpby(nrow, 1.0, x, 0.0, x, p);
        alp = 1.0;
        bta = 0.0;
        if (alp == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                p[i] = x[i] + bta * x[i];
        }
        else if (bta == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                p[i] = alp * x[i] + x[i];
        }
        else
        {
            for (int i = 0; i < nrow; i++)
                p[i] = alp * x[i] + bta * x[i];
        }
        // waxpy end

        // HPC_sparsemv(A, p, Ap);
        for (int i = 0; i < nrow; i++)
        {
            double sum = 0.0;

            for (int j = 0; j < A.nnz_in_row[i]; j++)
                sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
            Ap[i] = sum;
        }
        // HPC_sparsemv end

        // waxpby(nrow, 1.0, b, -1.0, Ap, r);
        alp = 1.0;
        bta = -1.0;
        if (alp == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                r[i] = b[i] + bta * Ap[i];
        }
        else if (bta == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                r[i] = alp * b[i] + Ap[i];
        }
        else
        {
            for (int i = 0; i < nrow; i++)
                r[i] = alp * b[i] + bta * Ap[i];
        }
        // waxpby end

        // ddot(nrow, r, r, &rtrans);
        rtrans = 0.0;
        for (int i = 0; i < nrow; i++)
            rtrans += r[i] * r[i];
        // ddot end

        normr = sqrt(rtrans);

        if (rank == 0)
            cout << "Initial Residual = " << normr << endl;

        for (int k = 1; k < max_iter && normr > tolerance; k++)
        {
            if (k == 1)
            {
                // waxpby(nrow, 1.0, r, 0.0, r, p);
                alp = 1.0;
                bta = 0.0;
                if (alp == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = r[i] + bta * r[i];
                }
                else if (bta == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + r[i];
                }
                else
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + bta * r[i];
                }
                // waxpby end
            }
            else
            {
                oldrtrans = rtrans;

                // ddot (nrow, r, r, &rtrans); // 2*nrow ops
                rtrans = 0.0;
                for (int i = 0; i < nrow; i++)
                    rtrans += r[i] * r[i];
                // ddot end

                beta = rtrans / oldrtrans;

                // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
                alp = 1.0;
                bta = beta;
                if (alp == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = r[i] + bta * p[i];
                }
                else if (bta == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + p[i];
                }
                else
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + bta * p[i];
                }
                // waxpby end
            }

            normr = sqrt(rtrans);
            if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
                cout << "Iteration = " << k << "   Residual = " << normr << endl;

            // HPC_sparsemv(A, p, Ap); // 2*nnz ops
            for (int i = 0; i < nrow; i++)
            {
                double sum = 0.0;

                for (int j = 0; j < A.nnz_in_row[i]; j++)
                    sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
                Ap[i] = sum;
            }
            // HPC_sparsemv end
            for (int i = 0; i < nrow; i++)
                alpha += p[i] * Ap[i];
            // ddot end

            alpha = rtrans / alpha;

            // waxpby(nrow, 1.0, x, alpha, p, x);   // 2*nrow ops
            alp = 1.0;
            bta = alpha;
            if (alp == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = x[i] + bta * p[i];
            }
            else if (bta == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = alp * x[i] + p[i];
            }
            else
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = alp * x[i] + bta * p[i];
            }
            // waxpby end

            // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
            alp = 1.0;
            bta = -alpha;
            if (alp == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = r[i] + bta * Ap[i];
            }
            else if (bta == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = alp * r[i] + Ap[i];
            }
            else
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = alp * r[i] + bta * Ap[i];
            }
            // waxpby end

            niters = k;
        }

        // Compute difference between known exact solution and computed solution
        // All processors are needed here.

        float residual = 0;

        // compute_residual(A.local_nrow, x, xexact, &residual);
        for (int i = 0; i < nrow; i++)
        {
            float diff = fabs(x[i] - xexact[i]);
            if (diff > residual)
                residual = diff;
        }

        return residual;
    }
}

static void HPCCGLowerPrec(benchmark::State &state)
{
    double *x, *b, *xexact;
    double norm, d;
    int ierr = 0;
    int i, j;
    int ione = 1;

    int size = 1; // Serial case (not using MPI)
    int rank = 0;

    generate_matrix(nx, ny, nz, A, &x, &b, &xexact);

    bool dump_matrix = false;
    if (dump_matrix && size <= 4)
        dump_matlab_matrix(A);

    int niters = 0;
    float normr = 0.0;
    int max_iter = 100;
    float tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

    int nrow = A.local_nrow;
    int ncol = A.local_ncol;

    double *r = new double[nrow];
    float *p = new float[ncol]; // In parallel case, A is rectangular
    float *Ap = new float[nrow];
    float residual;

    cout_suppressor suppressor;
    for (auto _ : state)
    {
        residual =
            lower_prec::HPCCG_residual(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact);
    }

    delete[] p;
    delete[] Ap;
    delete[] r;

    if (rank == 0)
        cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
             << residual << ".\n"
             << endl;
}

namespace high_prec
{
    double HPCCG_residual(double *b, double *x, int max_iter, double tolerance,
                          int &niters, double &normr, double *r, double *p,
                          double *Ap, double *xexact)
    {
        // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
        int nrow = A.local_nrow;
        int ncol = A.local_ncol;

        normr = 0.0;
        double rtrans = 0.0;
        double oldrtrans = 0.0;
        double alp, bta, alpha, beta, diff;

        int rank = 0; // Serial case (not using MPI)

        int print_freq = max_iter / 10;
        if (print_freq > 50)
            print_freq = 50;
        if (print_freq < 1)
            print_freq = 1;

        // // p is of length ncols, copy x to p for sparse MV operation
        // waxpby(nrow, 1.0, x, 0.0, x, p);
        alp = 1.0;
        bta = 0.0;
        if (alp == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                p[i] = x[i] + bta * x[i];
        }
        else if (bta == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                p[i] = alp * x[i] + x[i];
        }
        else
        {
            for (int i = 0; i < nrow; i++)
                p[i] = alp * x[i] + bta * x[i];
        }
        // waxpy end

        // HPC_sparsemv(A, p, Ap);
        for (int i = 0; i < nrow; i++)
        {
            double sum = 0.0;

            for (int j = 0; j < A.nnz_in_row[i]; j++)
                sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
            Ap[i] = sum;
        }
        // HPC_sparsemv end

        // waxpby(nrow, 1.0, b, -1.0, Ap, r);
        alp = 1.0;
        bta = -1.0;
        if (alp == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                r[i] = b[i] + bta * Ap[i];
        }
        else if (bta == 1.0)
        {
            for (int i = 0; i < nrow; i++)
                r[i] = alp * b[i] + Ap[i];
        }
        else
        {
            for (int i = 0; i < nrow; i++)
                r[i] = alp * b[i] + bta * Ap[i];
        }
        // waxpby end

        // ddot(nrow, r, r, &rtrans);
        rtrans = 0.0;
        for (int i = 0; i < nrow; i++)
            rtrans += r[i] * r[i];
        // ddot end

        normr = sqrt(rtrans);

        if (rank == 0)
            cout << "Initial Residual = " << normr << endl;

        for (int k = 1; k < max_iter && normr > tolerance; k++)
        {
            if (k == 1)
            {
                // waxpby(nrow, 1.0, r, 0.0, r, p);
                alp = 1.0;
                bta = 0.0;
                if (alp == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = r[i] + bta * r[i];
                }
                else if (bta == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + r[i];
                }
                else
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + bta * r[i];
                }
                // waxpby end
            }
            else
            {
                oldrtrans = rtrans;

                // ddot (nrow, r, r, &rtrans); // 2*nrow ops
                rtrans = 0.0;
                for (int i = 0; i < nrow; i++)
                    rtrans += r[i] * r[i];
                // ddot end

                beta = rtrans / oldrtrans;

                // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
                alp = 1.0;
                bta = beta;
                if (alp == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = r[i] + bta * p[i];
                }
                else if (bta == 1.0)
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + p[i];
                }
                else
                {
                    for (int i = 0; i < nrow; i++)
                        p[i] = alp * r[i] + bta * p[i];
                }
                // waxpby end
            }

            normr = sqrt(rtrans);
            if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
                cout << "Iteration = " << k << "   Residual = " << normr << endl;

            // HPC_sparsemv(A, p, Ap); // 2*nnz ops
            for (int i = 0; i < nrow; i++)
            {
                double sum = 0.0;

                for (int j = 0; j < A.nnz_in_row[i]; j++)
                    sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
                Ap[i] = sum;
            }
            // HPC_sparsemv end

            alpha = 0.0;

            // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
            alpha = 0.0;
            for (int i = 0; i < nrow; i++)
                alpha += p[i] * Ap[i];
            // ddot end

            alpha = rtrans / alpha;

            // waxpby(nrow, 1.0, x, alpha, p, x);   // 2*nrow ops
            alp = 1.0;
            bta = alpha;
            if (alp == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = x[i] + bta * p[i];
            }
            else if (bta == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = alp * x[i] + p[i];
            }
            else
            {
                for (int i = 0; i < nrow; i++)
                    x[i] = alp * x[i] + bta * p[i];
            }
            // waxpby end

            // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
            alp = 1.0;
            bta = -alpha;
            if (alp == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = r[i] + bta * Ap[i];
            }
            else if (bta == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = alp * r[i] + Ap[i];
            }
            else
            {
                for (int i = 0; i < nrow; i++)
                    r[i] = alp * r[i] + bta * Ap[i];
            }
            // waxpby end

            niters = k;
        }

        // Compute difference between known exact solution and computed solution
        // All processors are needed here.

        double residual = 0;

        // compute_residual(A.local_nrow, x, xexact, &residual);
        for (int i = 0; i < nrow; i++)
        {
            diff = fabs(x[i] - xexact[i]);
            if (diff > residual)
                residual = diff;
        }

        return residual;
    }

} // high_prec

static void HPCCGHighPrec(benchmark::State &state)
{
    double *x, *b, *xexact;
    double norm, d;
    int ierr = 0;
    int i, j;
    int ione = 1;

    int size = 1; // Serial case (not using MPI)
    int rank = 0;

    generate_matrix(nx, ny, nz, A, &x, &b, &xexact);

    bool dump_matrix = false;
    if (dump_matrix && size <= 4)
        dump_matlab_matrix(A);

    int niters = 0;
    double normr = 0.0;
    int max_iter = 100;
    double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

    int nrow = A.local_nrow;
    int ncol = A.local_ncol;

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, A is rectangular
    double *Ap = new double[nrow];
    double residual;

    cout_suppressor suppressor;
    for (auto _ : state)
    {
        residual =
            high_prec::HPCCG_residual(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact);
    }

    delete[] p;
    delete[] Ap;
    delete[] r;

    if (rank == 0)
        cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
             << residual << ".\n"
             << endl;
}

BENCHMARK(HPCCGLowerPrec)->Unit(benchmark::kSecond)->Iterations(1);
BENCHMARK(HPCCGHighPrec)->Unit(benchmark::kSecond)->Iterations(1);

BENCHMARK_MAIN();