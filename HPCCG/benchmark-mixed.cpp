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
    void generate_matrix(HPC_Sparse_Matrix &A, double **x,
                         float **b, double **xexact)

    {
        int debug = 0;

        int size = 1; // Serial case (not using MPI)
        int rank = 0;

        // *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
        A.title = 0;

        // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
        bool use_7pt_stencil = false;

        int local_nrow = nx * ny * nz;   // This is the size of our subblock
        assert(local_nrow > 0);          // Must have something to work with
        int local_nnz = 27 * local_nrow; // Approximately 27 nonzeros per row (except
                                         // for boundary nodes)

        int total_nrow = local_nrow * size; // Total number of grid points in mesh
        long long total_nnz =
            27 * (long long)total_nrow; // Approximately 27 nonzeros per row (except
                                        // for boundary nodes)

        int start_row =
            local_nrow *
            rank; // Each processor gets a section of a chimney stack domain
        int stop_row = start_row + local_nrow - 1;

        // Allocate arrays that are of length local_nrow
        A.nnz_in_row = new int[local_nrow];
        A.ptr_to_vals_in_row = new double *[local_nrow];
        A.ptr_to_inds_in_row = new int *[local_nrow];
        A.ptr_to_diags = new double *[local_nrow];

        *x = new double[local_nrow];
        *b = new float[local_nrow];
        *xexact = new double[local_nrow];

        // Allocate arrays that are of length local_nnz
        A.list_of_vals = new double[local_nnz];
        A.list_of_inds = new int[local_nnz];

        double *curvalptr = A.list_of_vals;
        int *curindptr = A.list_of_inds;

        long long nnzglobal = 0;
        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int curlocalrow = iz * nx * ny + iy * nx + ix;
                    int currow = start_row + iz * nx * ny + iy * nx + ix;
                    int nnzrow = 0;
                    A.ptr_to_vals_in_row[curlocalrow] = curvalptr;
                    A.ptr_to_inds_in_row[curlocalrow] = curindptr;
                    for (int sz = -1; sz <= 1; sz++)
                    {
                        for (int sy = -1; sy <= 1; sy++)
                        {
                            for (int sx = -1; sx <= 1; sx++)
                            {
                                int curcol = currow + sz * nx * ny + sy * nx + sx;
                                //            Since we have a stack of nx by ny by nz domains ,
                                //            stacking in the z direction, we check to see if sx
                                //            and sy are reaching outside of the domain, while the
                                //            check for the curcol being valid is sufficient to
                                //            check the z values
                                if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) &&
                                    (iy + sy < ny) && (curcol >= 0 && curcol < total_nrow))
                                {
                                    if (!use_7pt_stencil ||
                                        (sz * sz + sy * sy + sx * sx <=
                                         1))
                                    { // This logic will skip over point that are not part
                                      // of a 7-pt stencil
                                        if (curcol == currow)
                                        {
                                            A.ptr_to_diags[curlocalrow] = curvalptr;
                                            *curvalptr++ = 27.0;
                                        }
                                        else
                                        {
                                            *curvalptr++ = -1.0;
                                        }
                                        *curindptr++ = curcol;
                                        nnzrow++;
                                    }
                                }
                            } // end sx loop
                        }     // end sy loop
                    }         // end sz loop
                    A.nnz_in_row[curlocalrow] = nnzrow;
                    nnzglobal += nnzrow;
                    (*x)[curlocalrow] = 0.0;
                    (*b)[curlocalrow] = 27.0 - ((double)(nnzrow - 1));
                    (*xexact)[curlocalrow] = 1.0;
                } // end ix loop
            }     // end iy loop
        }         // end iz loop
        if (debug)
            cout << "Process " << rank << " of " << size << " has " << local_nrow;

        if (debug)
            cout << " rows. Global rows " << start_row << " through " << stop_row
                 << endl;

        if (debug)
            cout << "Process " << rank << " of " << size << " has " << local_nnz
                 << " nonzeros." << endl;

        A.start_row = start_row;
        A.stop_row = stop_row;
        A.total_nrow = total_nrow;
        A.total_nnz = total_nnz;
        A.local_nrow = local_nrow;
        A.local_ncol = local_nrow;
        A.local_nnz = local_nnz;

        return;
    }

    double HPCCG(float *b, double *x, int max_iter, float tolerance,
                 int &niters, double &normr, double *r, double *p,
                 double *Ap, double *xexact)
    {
        // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
        int nrow = A.local_nrow;
        int ncol = A.local_ncol;

        normr = 0.0;
        double rtrans = 0.0;
        float rtransf = 0.0;
        double oldrtrans = 0.0;
        float oldrtransf = 0.0;
        float alp;
        double bta;
        float btaf;

        int rank = 0; // Serial case (not using MPI)

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
            for (int j = 0; j < A.nnz_in_row[i]; j++)
                Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
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

        int k = 1;

        for (; k < max_iter / 2; k++)
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
                oldrtrans = rtrans;
            }
            else
            {
                // ddot (nrow, r, r, &rtrans); // 2*nrow ops
                rtrans = 0.0;
                for (int i = 0; i < nrow; i++)
                    rtrans += r[i] * r[i];
                // ddot end

                double beta = rtrans / oldrtrans;

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
                oldrtrans = rtrans;
            }

            normr = sqrt(rtrans);

            // HPC_sparsemv(A, p, Ap); // 2*nnz ops
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < A.nnz_in_row[i]; j++)
                    Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
            }
            // HPC_sparsemv end

            double alpha = 0.0;

            // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
            alpha = 0.0;
            if (Ap == p)
                for (int i = 0; i < nrow; i++)
                    alpha += p[i] * p[i];
            else
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

        oldrtransf = oldrtrans;

        for (int)

        for (; k < max_iter; k++)
        {
            // ddot (nrow, r, r, &rtrans); // 2*nrow ops
            rtransf = 0.0;
            for (int i = 0; i < nrow; i++)
                rtransf += r[i] * r[i];
            // ddot end

            float beta = rtransf / oldrtransf;

            // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
            alp = 1.0;
            btaf = beta;
            if (alp == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    p[i] = r[i] + btaf * p[i];
            }
            else if (btaf == 1.0)
            {
                for (int i = 0; i < nrow; i++)
                    p[i] = alp * r[i] + p[i];
            }
            else
            {
                for (int i = 0; i < nrow; i++)
                    p[i] = alp * r[i] + btaf * p[i];
            }
            // waxpby end
            oldrtransf = rtransf;

            normr = sqrt(rtransf);

            // HPC_sparsemv(A, p, Ap); // 2*nnz ops
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < A.nnz_in_row[i]; j++)
                    Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
            }
            // HPC_sparsemv end

            float alpha = 0.0;

            // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
            alpha = 0.0;
            if (Ap == p)
                for (int i = 0; i < nrow; i++)
                    alpha += p[i] * p[i];
            else
                for (int i = 0; i < nrow; i++)
                    alpha += p[i] * Ap[i];
            // ddot end

            alpha = rtransf / alpha;

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

        // compute residual(A.local_nrow, x, xexact, &residual);
        for (int i = 0; i < nrow; i++)
        {
            float diff = fabs(x[i] - xexact[i]);
            if (diff > residual)
                residual = diff;
        }

        return residual;
    }

}

namespace high_prec
{
    double HPCCG(double const *b, double *x, const int max_iter, const double tolerance,
                 int &niters, double &normr, double *r, double *p,
                 double *Ap, double const *xexact)
    {
        // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
        int nrow = A.local_nrow;
        int ncol = A.local_ncol;

        normr = 0.0;
        double rtrans = 0.0;
        double oldrtrans = 0.0;
        double alp, bta;

        int rank = 0; // Serial case (not using MPI)

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
            for (int j = 0; j < A.nnz_in_row[i]; j++)
                Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
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

        for (int k = 1; k < max_iter; k++)
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

                double beta = rtrans / oldrtrans;

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

            // HPC_sparsemv(A, p, Ap); // 2*nnz ops
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < A.nnz_in_row[i]; j++)
                    Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
            }
            // HPC_sparsemv end

            double alpha = 0.0;

            // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
            alpha = 0.0;
            if (Ap == p)
                for (int i = 0; i < nrow; i++)
                    alpha += p[i] * p[i];
            else
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

        // compute residual(A.local_nrow, x, xexact, &residual);
        for (int i = 0; i < nrow; i++)
        {
            double diff = fabs(x[i] - xexact[i]);
            if (diff > residual)
                residual = diff;
        }

        return residual;
    }

} // high_prec

static void HPCCGLowerPrec(benchmark::State &state)
{
    double *x, *xexact;
    float *b;
    double norm, d;
    int ierr = 0;
    int i, j;
    int ione = 1;

    int size = 1; // Serial case (not using MPI)
    int rank = 0;

    lower_prec::generate_matrix(A, &x, &b, &xexact);

    bool dump_matrix = false;
    if (dump_matrix && size <= 4)
        dump_matlab_matrix(A);

    int nrow = A.local_nrow;
    int ncol = A.local_ncol;

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, A is rectangular
    double *Ap = new double[nrow];
    double *rf = new double[nrow];
    double *pf = new double[ncol]; // In parallel case, A is rectangular
    float *Apf = new float[nrow];
    float residual;

    int niters = 0;
    double normr = 0.0;
    int max_iter = 100;
    double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

    cout_suppressor suppressor;
    for (auto _ : state)
    {
        residual =
            lower_prec::HPCCG(b, x, max_iter, tolerance, niters, normr, r, p, Ap, rf, pf, Apf, xexact);
    
        for (int i = 0; i < nrow; i++) {
            x[i] = 0;
            rf[i] = 0;
            Apf[i] = 0;
        }
        for (int i = 0; i < ncol; i++) {
          pf[i] = 0;
        }
        niters = 0;
        normr = 0;
    }

    delete[] p;
    delete[] Ap;
    delete[] xexact;
    delete[] b;
    delete[] x;

    if (rank == 0)
        cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
             << residual << ".\n"
             << endl;
}

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

    int nrow = A.local_nrow;
    int ncol = A.local_ncol;

    double *r = new double[nrow];
    double *p = new double[ncol]; // In parallel case, A is rectangular
    double *Ap = new double[nrow];
    double residual;

    int niters = 0;
    double normr = 0.0;
    int max_iter = 100;
    double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

    cout_suppressor suppressor;
    for (auto _ : state)
    {
        residual =
            high_prec::HPCCG(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact);
        
        for (int i = 0; i < nrow; i++) {
            x[i] = 0;
            r[i] = 0;
            Ap[i] = 0;
        }
        for (int i = 0; i < ncol; i++) {
          p[i] = 0;
        }
        niters = 0;
        normr = 0;
    }

    delete[] p;
    delete[] Ap;
    delete[] xexact;
    delete[] b;
    delete[] x;
    delete[] r;

    if (rank == 0)
        cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
             << residual << ".\n"
             << endl;
}

BENCHMARK(HPCCGLowerPrec)->Unit(benchmark::kSecond);
BENCHMARK(HPCCGHighPrec)->Unit(benchmark::kSecond);

BENCHMARK_MAIN();
