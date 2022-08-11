#include <cmath>

#include "adapt.h"
#include "adapt-impl.cpp"

namespace adapt
{

    AD_real BlkSchlsEqEuroNoDiv(AD_real sptprice,
                                AD_real strike, AD_real rate, AD_real volatility,
                                AD_real time, int otype)
    {
        AD_INDEPENDENT(sptprice, "sptprice");
        AD_INDEPENDENT(strike, "strike");
        AD_INDEPENDENT(rate, "rate");
        AD_INDEPENDENT(volatility, "volatility");
        AD_INDEPENDENT(time, "time");

        AD_real OptionPrice;

        // local private working variables for the calculation
        AD_real xStockPrice;
        AD_real xStrikePrice;
        AD_real xRiskFreeRate;
        AD_real xVolatility;
        AD_real xTime;
        AD_real xSqrtTime;

        AD_real logValues;
        AD_real xLogTerm;
        AD_real xD1;
        AD_real xD2;
        AD_real xPowerTerm;
        AD_real xDen;
        AD_real d1;
        AD_real d2;
        AD_real FutureValueX;
        AD_real NofXd1;
        AD_real NofXd2;
        AD_real NegNofXd1;
        AD_real NegNofXd2;

        AD_real InputX;

        int sign;

        AD_real OutputX;
        AD_real xInput;
        AD_real xNPrimeofX;
        AD_real expValues;
        AD_real xK2;
        AD_real xK2_2, xK2_3;
        AD_real xK2_4, xK2_5;
        AD_real xLocal, xLocal_1;
        AD_real xLocal_2, xLocal_3;
        AD_real inv_sqrt_2xPI = 0.39894228040143270286;

        xStockPrice = sptprice;
        xStrikePrice = strike;
        xRiskFreeRate = rate;
        xVolatility = volatility;

        AD_INTERMEDIATE(xStockPrice, "xStockPrice");
        AD_INTERMEDIATE(xStrikePrice, "xStrikePrice");
        AD_INTERMEDIATE(xRiskFreeRate, "xRiskFreeRate");
        AD_INTERMEDIATE(xVolatility, "xVolatility");

        xTime = time;
        AD_INTERMEDIATE(xTime, "xTime");

        xSqrtTime = sqrt(xTime);
        AD_INTERMEDIATE(xSqrtTime, "xSqrtTime");

        logValues = log(sptprice / strike);
        AD_INTERMEDIATE(logValues, "logValues");

        xLogTerm = logValues;
        AD_INTERMEDIATE(xLogTerm, "xLogTerm");

        xPowerTerm = xVolatility * xVolatility;
        AD_INTERMEDIATE(xPowerTerm, "xPowerTerm");
        xPowerTerm = xPowerTerm * 0.5;
        AD_INTERMEDIATE(xPowerTerm, "xPowerTerm");

        xD1 = xRiskFreeRate + xPowerTerm;
        AD_INTERMEDIATE(xD1, "xD1");
        xD1 = xD1 * xTime;
        AD_INTERMEDIATE(xD1, "xD1");
        xD1 = xD1 + xLogTerm;
        AD_INTERMEDIATE(xD1, "xD1");

        xDen = xVolatility * xSqrtTime;
        AD_INTERMEDIATE(xDen, "xDen");
        xD1 = xD1 / xDen;
        AD_INTERMEDIATE(xD1, "xD1");
        xD2 = xD1 - xDen;
        AD_INTERMEDIATE(xD2, "xD2");

        // ---------------------------- NofXd1 = CNDF(d1) --------------------------
        InputX = xD1;

        // Check for negative value of InputX
        if (InputX < 0.0)
        {
            InputX = -InputX;
            sign = 1;
        }
        else
            sign = 0;

        AD_INTERMEDIATE(InputX, "InputX");

        xInput = InputX;
        AD_INTERMEDIATE(xInput, "xInput");

        // Compute NPrimeX term common to both four & six decimal accuracy calcs
        expValues = exp(-0.5f * InputX * InputX);
        AD_INTERMEDIATE(expValues, "expValues");
        xNPrimeofX = expValues;
        xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;
        AD_INTERMEDIATE(xNPrimeofX, "xNPrimeofX");

        xK2 = 0.2316419 * xInput;
        xK2 = 1.0 + xK2;
        xK2 = 1.0 / xK2;
        AD_INTERMEDIATE(xK2, "xK2");
        xK2_2 = xK2 * xK2;
        AD_INTERMEDIATE(xK2_2, "xK2_2");
        xK2_3 = xK2_2 * xK2;
        AD_INTERMEDIATE(xK2_3, "xK2_3");
        xK2_4 = xK2_3 * xK2;
        AD_INTERMEDIATE(xK2_4, "xK2_4");
        xK2_5 = xK2_4 * xK2;
        AD_INTERMEDIATE(xK2_5, "xK2_5");

        xLocal_1 = xK2 * 0.319381530;
        AD_INTERMEDIATE(xLocal_1, "xLocal_2");
        xLocal_2 = xK2_2 * (-0.356563782);
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");
        xLocal_3 = xK2_3 * 1.781477937;
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        xLocal_3 = xK2_4 * (-1.821255978);
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        xLocal_3 = xK2_5 * 1.330274429;
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");

        xLocal_1 = xLocal_2 + xLocal_1;
        AD_INTERMEDIATE(xLocal_1, "xLocal_1");
        xLocal = xLocal_1 * xNPrimeofX;
        AD_INTERMEDIATE(xLocal, "xLocal");
        xLocal = 1.0 - xLocal;
        AD_INTERMEDIATE(xLocal, "xLocal");

        OutputX = xLocal;

        if (sign)
        {
            OutputX = 1.0 - OutputX;
        }
        AD_INTERMEDIATE(OutputX, "OutputX");

        NofXd1 = OutputX;
        AD_INTERMEDIATE(NofXd1, "NofXd1");
        // -------------------------- NofXd1 = CNDF(d1) END ------------------------

        // ---------------------------- NofXd2 = CNDF(d2) --------------------------
        InputX = xD2;

        // Check for negative value of InputX
        if (InputX < 0.0)
        {
            InputX = -InputX;
            sign = 1;
        }
        else
            sign = 0;

        AD_INTERMEDIATE(InputX, "InputX");

        xInput = InputX;
        AD_INTERMEDIATE(xInput, "xInput");

        // Compute NPrimeX term common to both four & six decimal accuracy calcs
        expValues = exp(-0.5f * InputX * InputX);
        AD_INTERMEDIATE(expValues, "expValues");
        xNPrimeofX = expValues;
        xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;
        AD_INTERMEDIATE(xNPrimeofX, "xNPrimeofX");

        xK2 = 0.2316419 * xInput;
        xK2 = 1.0 + xK2;
        xK2 = 1.0 / xK2;
        AD_INTERMEDIATE(xK2, "xK2");
        xK2_2 = xK2 * xK2;
        AD_INTERMEDIATE(xK2_2, "xK2_2");
        xK2_3 = xK2_2 * xK2;
        AD_INTERMEDIATE(xK2_3, "xK2_3");
        xK2_4 = xK2_3 * xK2;
        AD_INTERMEDIATE(xK2_4, "xK2_4");
        xK2_5 = xK2_4 * xK2;
        AD_INTERMEDIATE(xK2_5, "xK2_5");

        xLocal_1 = xK2 * 0.319381530;
        AD_INTERMEDIATE(xLocal_1, "xLocal_1");
        xLocal_2 = xK2_2 * (-0.356563782);
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");
        xLocal_3 = xK2_3 * 1.781477937;
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");
        xLocal_3 = xK2_4 * (-1.821255978);
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");
        xLocal_3 = xK2_5 * 1.330274429;
        AD_INTERMEDIATE(xLocal_3, "xLocal_3");
        xLocal_2 = xLocal_2 + xLocal_3;
        AD_INTERMEDIATE(xLocal_2, "xLocal_2");

        xLocal_1 = xLocal_2 + xLocal_1;
        AD_INTERMEDIATE(xLocal_1, "xLocal_1");
        xLocal = xLocal_1 * xNPrimeofX;
        AD_INTERMEDIATE(xLocal, "xLocal");
        xLocal = 1.0 - xLocal;
        AD_INTERMEDIATE(xLocal, "xLocal");

        OutputX = xLocal;

        if (sign)
        {
            OutputX = 1.0 - OutputX;
        }
        AD_INTERMEDIATE(OutputX, "OutputX");

        NofXd2 = OutputX;
        AD_INTERMEDIATE(NofXd2, "NofXd2");
        // -------------------------- NofXd2 = CNDF(d2) END ------------------------

        FutureValueX = strike * (exp(-(rate) * (time)));
        AD_INTERMEDIATE(FutureValueX, "FutureValueX");
        if (otype == 0)
        {
            OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
        }
        else
        {
            NegNofXd1 = (1.0 - NofXd1);
            AD_INTERMEDIATE(NegNofXd1, "NegNofXd1");
            NegNofXd2 = (1.0 - NofXd2);
            AD_INTERMEDIATE(NegNofXd2, "NegNofXd2");
            OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
        }

        AD_DEPENDENT(OptionPrice, "OptionPrice", 0.0);

        return OptionPrice;
    }

}
