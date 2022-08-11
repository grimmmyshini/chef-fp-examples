#include <cstdio>
#include <cstring>

#include "benchmark/benchmark.h"

#include "blackscholes.hpp"
#include "blackscholes-adapt.hpp"

#include "adapt.h"

#define INPUT_FILE_NAME "data/adaptvsclad.txt"

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

static void ErrorEstimateBlkSolAdapt(benchmark::State &state)
{
    FILE *file;
    int i;
    int loopnum;
    AD_real *buffer;
    int *buffer2;
    int rv;

    OptionData *data;
    AD_real *prices;
    int numOptions;

    int *otype;
    AD_real *sptprice;
    AD_real *strike;
    AD_real *rate;
    AD_real *volatility;
    AD_real *otime;
    int numError = 0;
    int nThreads;

    nThreads = 1;
    const char *inputFile = INPUT_FILE_NAME;
    const char *outputFile = "prices.output";

    // Read input data from file
    file = fopen(inputFile, "r");
    if (file == NULL)
    {
        printf("ERROR: Unable to open file `%s'.\n", inputFile);
        exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if (rv != 1)
    {
        printf("ERROR: Unable to read from file `%s'.\n", inputFile);
        fclose(file);
        exit(1);
    }
    if (nThreads > numOptions)
    {
        printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
        nThreads = numOptions;
    }

    if (nThreads != 1)
    {
        printf("Error: <nthreads> must be 1 (serial version)\n");
        exit(1);
    }

    // alloc spaces for the option data
    data = new OptionData[numOptions];
    prices = new AD_real[numOptions];
    for (loopnum = 0; loopnum < numOptions; ++loopnum)
    {
        rv = fscanf(file, INPUT_LINE_FORMAT,
                    &data[loopnum].s,
                    &data[loopnum].strike,
                    &data[loopnum].r,
                    &data[loopnum].divq,
                    &data[loopnum].v,
                    &data[loopnum].t,
                    &data[loopnum].OptionType,
                    &data[loopnum].divs,
                    &data[loopnum].DGrefval);

        if (rv != 9)
        {
            printf("ERROR: Unable to read from file `%s'.\n", inputFile);
            fclose(file);
            exit(1);
        }
    }
    rv = fclose(file);
    if (rv != 0)
    {
        printf("ERROR: Unable to close file `%s'.\n", inputFile);
        exit(1);
    }

    const int LINESIZE = 64;
    const int PAD_FPTYPE = LINESIZE * sizeof(fptype);
    const int PAD_INT = LINESIZE * sizeof(int);

    buffer = new AD_real[5 * numOptions + LINESIZE];
    sptprice = (AD_real *)(((unsigned long long)buffer + PAD_FPTYPE) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = new int[numOptions + LINESIZE];
    otype = (int *)(((unsigned long long)buffer2 + PAD_INT) & ~(LINESIZE - 1));

    for (i = 0; i < numOptions; i++)
    {
        otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i] = data[i].s;
        strike[i] = data[i].strike;
        rate[i] = data[i].r;
        volatility[i] = data[i].v;
        otime[i] = data[i].t;
    }

    // serial version
    // ---------------------------- bs_thread() ------------------------------
    int k, j;
    AD_real price;
    fptype priceDelta;
    int start = 0 * (numOptions / nThreads);
    int end = start + (numOptions / nThreads);

    cout_suppressor suppressor;

    for (auto _ : state)
    {
        AD_begin();
        AD_enable_source_aggregation();

        for (j = 0; j < NUM_RUNS; j++)
        {
            for (k = start; k < end; k++)
            {
                /* Calling main function to calculate option value based on
                 * Black & Scholes's equation.
                 */
                price = adapt::BlkSchlsEqEuroNoDiv(sptprice[k], strike[k],
                                                   rate[k], volatility[k], otime[k],
                                                   otype[k]);
                prices[k] = price;

#ifdef ERR_CHK
                priceDelta = data[k].DGrefval - price;
                if (fabs(priceDelta) >= 1e-4)
                {
                    printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                           k, price, data[k].DGrefval, priceDelta);
                    numError++;
                }
#endif
            }
        }

        AD_report();
    }
    // -------------------------- bs_thread END -----------------------------


    delete[] buffer2;
    delete[] buffer;
    delete[] data;
    delete[] prices;
}

BENCHMARK(ErrorEstimateBlkSolAdapt)->Unit(benchmark::kSecond);

BENCHMARK_MAIN();