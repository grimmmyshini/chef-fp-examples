#include <cmath>
#include "fastonebigheader.h"

typedef double fptype;
#define INPUT_LINE_FORMAT "%lf %lf %lf %lf %lf %lf %c %lf %lf"

#define NUM_RUNS 100

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
fptype ApproxBlkSchlsEqEuroNoDiv(fptype sptprice,
                           fptype strike, fptype rate, fptype volatility,
                           fptype time, int otype)
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;

    fptype InputX;

    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;
    fptype inv_sqrt_2xPI = 0.39894228040143270286;

    fptype clad_log_logValues_;
    fptype clad_sqr_xSqrtTime_;
    fptype clad_exp_expValues_;
    fptype clad_exp_FutureValueX_;

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    clad_sqr_xSqrtTime_ = xTime;
    xSqrtTime = fastpow(clad_sqr_xSqrtTime_, 0.5);

    clad_log_logValues_ = sptprice / strike;
    logValues = fastlog(clad_log_logValues_);

    xLogTerm = logValues;

    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 - xDen;

    d1 = xD1;
    d2 = xD2;

    // ---------------------------- NofXd1 = CNDF(d1) --------------------------
    InputX = d1;

    // Check for negative value of InputX
    if (InputX < 0.0)
    {
        InputX = -InputX;
        sign = 1;
    }
    else
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    clad_exp_expValues_ = -0.5f * InputX * InputX;
    #ifdef FASTEXP
    expValues = fastexp(clad_exp_expValues_);
    #else
    expValues = exp(clad_exp_expValues_);
    #endif 
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal = xLocal_1 * xNPrimeofX;
    xLocal = 1.0 - xLocal;

    OutputX = xLocal;

    if (sign)
    {
        OutputX = 1.0 - OutputX;
    }

    NofXd1 = OutputX;
    // -------------------------- NofXd1 = CNDF(d1) END ------------------------

    // ---------------------------- NofXd2 = CNDF(d2) --------------------------
    InputX = d2;

    // Check for negative value of InputX
    if (InputX < 0.0)
    {
        InputX = -InputX;
        sign = 1;
    }
    else
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    clad_exp_expValues_ = -0.5f * InputX * InputX;
    #ifdef FASTEXP
    expValues = fastexp(clad_exp_expValues_);
    #else
    expValues = exp(clad_exp_expValues_);
    #endif 
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal = xLocal_1 * xNPrimeofX;
    xLocal = 1.0 - xLocal;

    OutputX = xLocal;

    if (sign)
    {
        OutputX = 1.0 - OutputX;
    }

    NofXd2 = OutputX;
    // -------------------------- NofXd2 = CNDF(d2) END ------------------------

    clad_exp_FutureValueX_ = -(rate) * (time);
    #ifdef FASTEXP
    FutureValueX = strike * (fastexp(clad_exp_FutureValueX_));
    #else
    FutureValueX = strike * (exp(clad_exp_FutureValueX_));
    #endif
    if (otype == 0)
    {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    }
    else
    {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}
