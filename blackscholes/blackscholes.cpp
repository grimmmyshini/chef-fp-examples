// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice
// Hall, John C. Hull,

#include <cstdio>
#include <cstring>

#include "blackscholes.hpp"

// Precision to use for calculations

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    FILE *file;
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

    printf("PARSEC Benchmark Suite\n");
    fflush(NULL);

    if (argc != 2)
    {
        printf("Usage:\n\t%s <inputFile>\n", argv[0]);
        exit(1);
    }
    nThreads = 1;
    char *inputFile = argv[1];
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
    prices = new fptype[numOptions];
    for (loopnum = 0; loopnum < numOptions; ++loopnum)
    {
        rv = fscanf(file, "%lf %lf %lf %lf %lf %lf %c %lf %lf",
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

    printf("Num of Options: %d\n", numOptions);
    printf("Num of Runs: %d\n", NUM_RUNS);


    const int LINESIZE = 64;
    const int PAD = LINESIZE * sizeof(fptype);

    buffer = new fptype[5 * numOptions + LINESIZE];
    sptprice = (fptype *)(((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = new int[numOptions + LINESIZE];
    otype = (int *)(((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    for (i = 0; i < numOptions; i++)
    {
        otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i] = data[i].s;
        strike[i] = data[i].strike;
        rate[i] = data[i].r;
        volatility[i] = data[i].v;
        otime[i] = data[i].t;
    }

    printf("Size of data: %ld\n", numOptions * (sizeof(OptionData) + sizeof(int)));

    // serial version
    // ---------------------------- bs_thread() ------------------------------
    int k, j;
    fptype price;
    fptype priceDelta;
    int start = 0 * (numOptions / nThreads);
    int end = start + (numOptions / nThreads);

    for (j = 0; j < NUM_RUNS; j++)
    {
        for (k = start; k < end; k++)
        {
            /* Calling main function to calculate option value based on
             * Black & Scholes's equation.
             */
            price = BlkSchlsEqEuroNoDiv(sptprice[k], strike[k],
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
    // -------------------------- bs_thread END -----------------------------

    // Write prices to output file
    file = fopen(outputFile, "w");
    if (file == NULL)
    {
        printf("ERROR: Unable to open file `%s'.\n", outputFile);
        exit(1);
    }
    rv = fprintf(file, "%i\n", numOptions);
    if (rv < 0)
    {
        printf("ERROR: Unable to write to file `%s'.\n", outputFile);
        fclose(file);
        exit(1);
    }
    for (i = 0; i < numOptions; i++)
    {
        rv = fprintf(file, "%.18f\n", prices[i]);
        if (rv < 0)
        {
            printf("ERROR: Unable to write to file `%s'.\n", outputFile);
            fclose(file);
            exit(1);
        }
    }
    rv = fclose(file);
    if (rv != 0)
    {
        printf("ERROR: Unable to close file `%s'.\n", outputFile);
        exit(1);
    }

    delete[] buffer2;
    delete[] buffer;
    delete[] data;
    delete[] prices;

    return 0;
}
