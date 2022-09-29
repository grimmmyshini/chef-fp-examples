#include <cstdio>
#include <cstring>
#include <iostream>

#include "blackscholes-approx.hpp"
#include "blackscholes.hpp"

#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"

const char *input_files[] = {"data/100.txt", "data/1000.txt", "data/10000.txt", "data/100000.txt", "data/1000000.txt"};

bool print = true;
static void Blackscholes_MP(benchmark::State &state)
{
  FILE *file;
  int opt = state.range(0);
  int i;
  int loopnum;
  fptype *buffer;
  int *buffer2;
  int rv;

  OptionData *data;
  fptype *prices;
  int numOptions;

  int *otype;
  fptype *sptprice;
  fptype *strike;
  fptype *rate;
  fptype *volatility;
  fptype *otime;
  int numError = 0;
  int nThreads;

  nThreads = 1;
  const char *outputFile = "prices.output";

  const char *inputFile = input_files[state.range(1)];
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
    printf("WARNING: Not enough work, reducing number of threads to match "
           "number of options.\n");
    nThreads = numOptions;
  }

  if (nThreads != 1)
  {
    printf("Error: <nthreads> must be 1 (serial version)\n");
    exit(1);
  }

  // alloc spaces for the option data
  data = new OptionData[numOptions];
  prices = new fptype[numOptions];
  for (loopnum = 0; loopnum < numOptions; ++loopnum)
  {
    rv = fscanf(file, INPUT_LINE_FORMAT, &data[loopnum].s,
                &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq,
                &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType,
                &data[loopnum].divs, &data[loopnum].DGrefval);

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

  buffer = new fptype[5 * numOptions + LINESIZE];
  sptprice =
      (fptype *)(((unsigned long long)buffer + PAD_FPTYPE) & ~(LINESIZE - 1));
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

  if (print)
  {
    printf("Num of Options: %d\n", numOptions);
    printf("Num of Runs: %d\n", NUM_RUNS);
    printf("Size of data: %ld\n",
           numOptions * (sizeof(OptionData) + sizeof(int)));
  }

  double maxError = 0;
  int k, j;
  fptype price;
  fptype priceDelta;
  int start = 0 * (numOptions / nThreads);
  int end = start + (numOptions / nThreads);

  if (!opt)
  {
    for (auto _ : state)
    {
      for (j = 0; j < NUM_RUNS; j++)
      {
        for (k = start; k < end; k++)
        {
          prices[k] = BlkSchlsEqEuroNoDiv(sptprice[k], strike[k], rate[k],
                                          volatility[k], otime[k], otype[k]);
          benchmark::DoNotOptimize(prices[k]);
        }
      }
    }
  }
  else
  {
    for (auto _ : state)
    {
      for (j = 0; j < NUM_RUNS; j++)
      {
        for (k = start; k < end; k++)
        {
          prices[k] =
              ApproxBlkSchlsEqEuroNoDiv(sptprice[k], strike[k], rate[k],
                                        volatility[k], otime[k], otype[k]);
          benchmark::DoNotOptimize(prices[k]);
        }
      }
    }
  }
#ifdef MIXEDCOMP
  double totalError = 0;
  for (j = 0; j < NUM_RUNS; j++)
  {
    for (k = start; k < end; k++)
    {
      double priceApprox = ApproxBlkSchlsEqEuroNoDiv(
          sptprice[k], strike[k], rate[k], volatility[k], otime[k], otype[k]);

      double priceExact = BlkSchlsEqEuroNoDiv(
          sptprice[k], strike[k], rate[k], volatility[k], otime[k], otype[k]);
      double diff = std::fabs(priceApprox - priceExact);
      totalError += diff;
      maxError = std::max(maxError, diff);
    }
  }
  if (print)
  {
    printf("Maximum error between exact and approx is: %f\n", maxError);
    printf("Total error between exact and approx is: %f\n", totalError);
    printf("Average error between exact and approx is: %f\n", totalError / (NUM_RUNS * (end - start)));
  }
#endif

  print = false;
  delete[] buffer2;
  delete[] buffer;
  delete[] data;
  delete[] prices;
}

// Define our main
// Compile with -DFASTEXP=On to enable fast exp calls.
// use -DMIXEDCOMP=On to also run approx vs mixed comparision.

// 0 -> Exact blackscholes
// 1 -> Blackscholes with fastapprox
BENCHMARK(Blackscholes_MP)->ArgsProduct({benchmark::CreateDenseRange(0, 1, 1), benchmark::CreateDenseRange(0, 4, 1)});

int main(int argc, char **argv)
{
  ::benchmark::RegisterMemoryManager(mm.get());
  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::RegisterMemoryManager(nullptr);
}