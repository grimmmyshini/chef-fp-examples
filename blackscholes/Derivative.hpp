namespace clad
{
    void BlkSchlsEqEuroNoDiv_grad(float sptprice, float strike, float rate, float volatility, float time, int otype, clad::array_ref<float> _d_sptprice, clad::array_ref<float> _d_strike, clad::array_ref<float> _d_rate, clad::array_ref<float> _d_volatility, clad::array_ref<float> _d_time, clad::array_ref<int> _d_otype, double &_final_error)
    {
        float _d_OptionPrice = 0;
        double _delta_OptionPrice = 0;
        float _EERepl_OptionPrice0;
        float _d_xStockPrice = 0;
        double _delta_xStockPrice = 0;
        float _EERepl_xStockPrice0;
        float _d_xStrikePrice = 0;
        double _delta_xStrikePrice = 0;
        float _EERepl_xStrikePrice0;
        float _d_xRiskFreeRate = 0;
        double _delta_xRiskFreeRate = 0;
        float _EERepl_xRiskFreeRate0;
        float _d_xVolatility = 0;
        double _delta_xVolatility = 0;
        float _EERepl_xVolatility0;
        float _d_xTime = 0;
        double _delta_xTime = 0;
        float _EERepl_xTime0;
        float _d_xSqrtTime = 0;
        double _delta_xSqrtTime = 0;
        float _EERepl_xSqrtTime0;
        float _d_logValues = 0;
        double _delta_logValues = 0;
        float _EERepl_logValues0;
        float _d_xLogTerm = 0;
        double _delta_xLogTerm = 0;
        float _EERepl_xLogTerm0;
        float _d_xD1 = 0;
        double _delta_xD1 = 0;
        float _EERepl_xD10;
        float _d_xD2 = 0;
        double _delta_xD2 = 0;
        float _EERepl_xD20;
        float _d_xPowerTerm = 0;
        double _delta_xPowerTerm = 0;
        float _EERepl_xPowerTerm0;
        float _d_xDen = 0;
        double _delta_xDen = 0;
        float _EERepl_xDen0;
        float _d_d1 = 0;
        double _delta_d1 = 0;
        float _EERepl_d10;
        float _d_d2 = 0;
        double _delta_d2 = 0;
        float _EERepl_d20;
        float _d_FutureValueX = 0;
        double _delta_FutureValueX = 0;
        float _EERepl_FutureValueX0;
        float _d_NofXd1 = 0;
        double _delta_NofXd1 = 0;
        float _EERepl_NofXd10;
        float _d_NofXd2 = 0;
        double _delta_NofXd2 = 0;
        float _EERepl_NofXd20;
        float _d_NegNofXd1 = 0;
        double _delta_NegNofXd1 = 0;
        float _EERepl_NegNofXd10;
        float _d_NegNofXd2 = 0;
        double _delta_NegNofXd2 = 0;
        float _EERepl_NegNofXd20;
        float _d_InputX = 0;
        double _delta_InputX = 0;
        float _EERepl_InputX0;
        int _d_sign = 0;
        float _d_OutputX = 0;
        double _delta_OutputX = 0;
        float _EERepl_OutputX0;
        float _d_xInput = 0;
        double _delta_xInput = 0;
        float _EERepl_xInput0;
        float _d_xNPrimeofX = 0;
        double _delta_xNPrimeofX = 0;
        float _EERepl_xNPrimeofX0;
        float _d_expValues = 0;
        double _delta_expValues = 0;
        float _EERepl_expValues0;
        float _d_xK2 = 0;
        double _delta_xK2 = 0;
        float _EERepl_xK20;
        float _d_xK2_2 = 0, _d_xK2_3 = 0;
        double _delta_xK2_2 = 0;
        float _EERepl_xK2_20;
        double _delta_xK2_20 = 0;
        float _EERepl_xK2_21;
        float _d_xK2_4 = 0, _d_xK2_5 = 0;
        double _delta_xK2_4 = 0;
        float _EERepl_xK2_40;
        double _delta_xK2_40 = 0;
        float _EERepl_xK2_41;
        float _d_xLocal = 0, _d_xLocal_1 = 0;
        double _delta_xLocal = 0;
        float _EERepl_xLocal0;
        double _delta_xLocal0 = 0;
        float _EERepl_xLocal1;
        float _d_xLocal_2 = 0, _d_xLocal_3 = 0;
        double _delta_xLocal_2 = 0;
        float _EERepl_xLocal_20;
        double _delta_xLocal_20 = 0;
        float _EERepl_xLocal_21;
        float _d_inv_sqrt_2xPI = 0;
        double _delta_inv_sqrt_2xPI = 0;
        float _EERepl_inv_sqrt_2xPI0;
        float _t0;
        float _EERepl_xSqrtTime1;
        float _t1;
        float _t2;
        float _t3;
        float _EERepl_logValues1;
        float _t4;
        float _t5;
        float _EERepl_xPowerTerm1;
        float _EERepl_xPowerTerm2;
        float _EERepl_xD11;
        float _t6;
        float _t7;
        float _EERepl_xD12;
        float _EERepl_xD13;
        float _t8;
        float _t9;
        float _EERepl_xDen1;
        float _t10;
        float _t11;
        float _EERepl_xD14;
        float _EERepl_xD21;
        bool _cond0;
        float _EERepl_InputX1;
        float _t12;
        float _t13;
        float _t14;
        float _t15;
        float _EERepl_expValues1;
        float _t16;
        float _t17;
        float _EERepl_xNPrimeofX1;
        double _t18;
        float _EERepl_xK21;
        float _EERepl_xK22;
        double _t19;
        float _EERepl_xK23;
        float _t20;
        float _t21;
        float _EERepl_xK2_22;
        double _delta_xK2_3 = 0;
        float _t22;
        float _t23;
        float _EERepl_xK2_31;
        float _t24;
        float _t25;
        float _EERepl_xK2_42;
        double _delta_xK2_5 = 0;
        float _t26;
        float _t27;
        float _EERepl_xK2_51;
        double _delta_xLocal_1 = 0;
        float _EERepl_xLocal_11;
        float _EERepl_xLocal_22;
        double _delta_xLocal_3 = 0;
        float _EERepl_xLocal_31;
        float _EERepl_xLocal_23;
        float _EERepl_xLocal_32;
        float _EERepl_xLocal_24;
        float _EERepl_xLocal_33;
        float _EERepl_xLocal_25;
        float _EERepl_xLocal_12;
        float _t28;
        float _t29;
        float _EERepl_xLocal2;
        float _EERepl_xLocal3;
        int _cond1;
        float _EERepl_OutputX1;
        bool _cond2;
        float _EERepl_InputX2;
        float _t30;
        float _t31;
        float _t32;
        float _t33;
        float _EERepl_expValues2;
        float _t34;
        float _t35;
        float _EERepl_xNPrimeofX2;
        double _t36;
        float _EERepl_xK24;
        float _EERepl_xK25;
        double _t37;
        float _EERepl_xK26;
        float _t38;
        float _t39;
        float _EERepl_xK2_23;
        float _t40;
        float _t41;
        float _EERepl_xK2_32;
        float _t42;
        float _t43;
        float _EERepl_xK2_43;
        float _t44;
        float _t45;
        float _EERepl_xK2_52;
        float _EERepl_xLocal_13;
        float _EERepl_xLocal_26;
        float _EERepl_xLocal_34;
        float _EERepl_xLocal_27;
        float _EERepl_xLocal_35;
        float _EERepl_xLocal_28;
        float _EERepl_xLocal_36;
        float _EERepl_xLocal_29;
        float _EERepl_xLocal_14;
        float _t46;
        float _t47;
        float _EERepl_xLocal4;
        float _EERepl_xLocal5;
        int _cond3;
        float _EERepl_OutputX2;
        double _t48;
        float _t49;
        float _t50;
        float _t51;
        float _t52;
        float _EERepl_FutureValueX1;
        bool _cond4;
        float _t53;
        float _t54;
        float _t55;
        float _t56;
        float _EERepl_OptionPrice1;
        float _EERepl_NegNofXd11;
        float _EERepl_NegNofXd21;
        float _t57;
        float _t58;
        float _t59;
        float _t60;
        float _EERepl_OptionPrice2;
        float OptionPrice;
        _EERepl_OptionPrice0 = OptionPrice;
        float xStockPrice;
        _EERepl_xStockPrice0 = xStockPrice;
        float xStrikePrice;
        _EERepl_xStrikePrice0 = xStrikePrice;
        float xRiskFreeRate;
        _EERepl_xRiskFreeRate0 = xRiskFreeRate;
        float xVolatility;
        _EERepl_xVolatility0 = xVolatility;
        float xTime;
        _EERepl_xTime0 = xTime;
        float xSqrtTime;
        _EERepl_xSqrtTime0 = xSqrtTime;
        float logValues;
        _EERepl_logValues0 = logValues;
        float xLogTerm;
        _EERepl_xLogTerm0 = xLogTerm;
        float xD1;
        _EERepl_xD10 = xD1;
        float xD2;
        _EERepl_xD20 = xD2;
        float xPowerTerm;
        _EERepl_xPowerTerm0 = xPowerTerm;
        float xDen;
        _EERepl_xDen0 = xDen;
        float d1;
        _EERepl_d10 = d1;
        float d2;
        _EERepl_d20 = d2;
        float FutureValueX;
        _EERepl_FutureValueX0 = FutureValueX;
        float NofXd1;
        _EERepl_NofXd10 = NofXd1;
        float NofXd2;
        _EERepl_NofXd20 = NofXd2;
        float NegNofXd1;
        _EERepl_NegNofXd10 = NegNofXd1;
        float NegNofXd2;
        _EERepl_NegNofXd20 = NegNofXd2;
        float InputX;
        _EERepl_InputX0 = InputX;
        int sign;
        float OutputX;
        _EERepl_OutputX0 = OutputX;
        float xInput;
        _EERepl_xInput0 = xInput;
        float xNPrimeofX;
        _EERepl_xNPrimeofX0 = xNPrimeofX;
        float expValues;
        _EERepl_expValues0 = expValues;
        float xK2;
        _EERepl_xK20 = xK2;
        float xK2_2, xK2_3;
        _EERepl_xK2_21 = xK2_2;
        _EERepl_xK2_20 = xK2_2;
        float xK2_4, xK2_5;
        _EERepl_xK2_41 = xK2_4;
        _EERepl_xK2_40 = xK2_4;
        float xLocal, xLocal_1;
        _EERepl_xLocal1 = xLocal;
        _EERepl_xLocal0 = xLocal;
        float xLocal_2, xLocal_3;
        _EERepl_xLocal_21 = xLocal_2;
        _EERepl_xLocal_20 = xLocal_2;
        float inv_sqrt_2xPI = 0.3989422804014327;
        _EERepl_inv_sqrt_2xPI0 = inv_sqrt_2xPI;
        xStockPrice = sptprice;
        xStrikePrice = strike;
        xRiskFreeRate = rate;
        xVolatility = volatility;
        xTime = time;
        _t0 = xTime;
        xSqrtTime = sqrt(_t0);
        _EERepl_xSqrtTime1 = xSqrtTime;
        _t2 = sptprice;
        _t1 = strike;
        _t3 = _t2 / _t1;
        logValues = log(_t3);
        _EERepl_logValues1 = logValues;
        xLogTerm = logValues;
        _t5 = xVolatility;
        _t4 = xVolatility;
        xPowerTerm = _t5 * _t4;
        _EERepl_xPowerTerm1 = xPowerTerm;
        xPowerTerm = xPowerTerm * 0.5;
        _EERepl_xPowerTerm2 = xPowerTerm;
        xD1 = xRiskFreeRate + xPowerTerm;
        _EERepl_xD11 = xD1;
        _t7 = xD1;
        _t6 = xTime;
        xD1 = _t7 * _t6;
        _EERepl_xD12 = xD1;
        xD1 = xD1 + xLogTerm;
        _EERepl_xD13 = xD1;
        _t9 = xVolatility;
        _t8 = xSqrtTime;
        xDen = _t9 * _t8;
        _EERepl_xDen1 = xDen;
        _t11 = xD1;
        _t10 = xDen;
        xD1 = _t11 / _t10;
        _EERepl_xD14 = xD1;
        xD2 = xD1 - xDen;
        _EERepl_xD21 = xD2;
        d1 = xD1;
        d2 = xD2;
        InputX = d1;
        _cond0 = InputX < 0.;
        if (_cond0)
        {
            InputX = -InputX;
            _EERepl_InputX1 = InputX;
            sign = 1;
        }
        else
            sign = 0;
        xInput = InputX;
        _t13 = InputX;
        _t14 = -0.5F * _t13;
        _t12 = InputX;
        _t15 = _t14 * _t12;
        expValues = exp(_t15);
        _EERepl_expValues1 = expValues;
        xNPrimeofX = expValues;
        _t17 = xNPrimeofX;
        _t16 = inv_sqrt_2xPI;
        xNPrimeofX = _t17 * _t16;
        _EERepl_xNPrimeofX1 = xNPrimeofX;
        _t18 = xInput;
        xK2 = 0.23164190000000001 * _t18;
        _EERepl_xK21 = xK2;
        xK2 = 1. + xK2;
        _EERepl_xK22 = xK2;
        _t19 = xK2;
        xK2 = 1. / _t19;
        _EERepl_xK23 = xK2;
        _t21 = xK2;
        _t20 = xK2;
        xK2_2 = _t21 * _t20;
        _EERepl_xK2_22 = xK2_2;
        _t23 = xK2_2;
        _t22 = xK2;
        xK2_3 = _t23 * _t22;
        _EERepl_xK2_31 = xK2_3;
        _t25 = xK2_3;
        _t24 = xK2;
        xK2_4 = _t25 * _t24;
        _EERepl_xK2_42 = xK2_4;
        _t27 = xK2_4;
        _t26 = xK2;
        xK2_5 = _t27 * _t26;
        _EERepl_xK2_51 = xK2_5;
        xLocal_1 = xK2 * 0.31938153000000002;
        _EERepl_xLocal_11 = xLocal_1;
        xLocal_2 = xK2_2 * (-0.356563782);
        _EERepl_xLocal_22 = xLocal_2;
        xLocal_3 = xK2_3 * 1.781477937;
        _EERepl_xLocal_31 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_23 = xLocal_2;
        xLocal_3 = xK2_4 * (-1.8212559779999999);
        _EERepl_xLocal_32 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_24 = xLocal_2;
        xLocal_3 = xK2_5 * 1.3302744289999999;
        _EERepl_xLocal_33 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_25 = xLocal_2;
        xLocal_1 = xLocal_2 + xLocal_1;
        _EERepl_xLocal_12 = xLocal_1;
        _t29 = xLocal_1;
        _t28 = xNPrimeofX;
        xLocal = _t29 * _t28;
        _EERepl_xLocal2 = xLocal;
        xLocal = 1. - xLocal;
        _EERepl_xLocal3 = xLocal;
        OutputX = xLocal;
        _cond1 = sign;
        if (_cond1)
        {
            OutputX = 1. - OutputX;
            _EERepl_OutputX1 = OutputX;
        }
        NofXd1 = OutputX;
        InputX = d2;
        _cond2 = InputX < 0.;
        if (_cond2)
        {
            InputX = -InputX;
            _EERepl_InputX2 = InputX;
            sign = 1;
        }
        else
            sign = 0;
        xInput = InputX;
        _t31 = InputX;
        _t32 = -0.5F * _t31;
        _t30 = InputX;
        _t33 = _t32 * _t30;
        expValues = exp(_t33);
        _EERepl_expValues2 = expValues;
        xNPrimeofX = expValues;
        _t35 = xNPrimeofX;
        _t34 = inv_sqrt_2xPI;
        xNPrimeofX = _t35 * _t34;
        _EERepl_xNPrimeofX2 = xNPrimeofX;
        _t36 = xInput;
        xK2 = 0.23164190000000001 * _t36;
        _EERepl_xK24 = xK2;
        xK2 = 1. + xK2;
        _EERepl_xK25 = xK2;
        _t37 = xK2;
        xK2 = 1. / _t37;
        _EERepl_xK26 = xK2;
        _t39 = xK2;
        _t38 = xK2;
        xK2_2 = _t39 * _t38;
        _EERepl_xK2_23 = xK2_2;
        _t41 = xK2_2;
        _t40 = xK2;
        xK2_3 = _t41 * _t40;
        _EERepl_xK2_32 = xK2_3;
        _t43 = xK2_3;
        _t42 = xK2;
        xK2_4 = _t43 * _t42;
        _EERepl_xK2_43 = xK2_4;
        _t45 = xK2_4;
        _t44 = xK2;
        xK2_5 = _t45 * _t44;
        _EERepl_xK2_52 = xK2_5;
        xLocal_1 = xK2 * 0.31938153000000002;
        _EERepl_xLocal_13 = xLocal_1;
        xLocal_2 = xK2_2 * (-0.356563782);
        _EERepl_xLocal_26 = xLocal_2;
        xLocal_3 = xK2_3 * 1.781477937;
        _EERepl_xLocal_34 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_27 = xLocal_2;
        xLocal_3 = xK2_4 * (-1.8212559779999999);
        _EERepl_xLocal_35 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_28 = xLocal_2;
        xLocal_3 = xK2_5 * 1.3302744289999999;
        _EERepl_xLocal_36 = xLocal_3;
        xLocal_2 = xLocal_2 + xLocal_3;
        _EERepl_xLocal_29 = xLocal_2;
        xLocal_1 = xLocal_2 + xLocal_1;
        _EERepl_xLocal_14 = xLocal_1;
        _t47 = xLocal_1;
        _t46 = xNPrimeofX;
        xLocal = _t47 * _t46;
        _EERepl_xLocal4 = xLocal;
        xLocal = 1. - xLocal;
        _EERepl_xLocal5 = xLocal;
        OutputX = xLocal;
        _cond3 = sign;
        if (_cond3)
        {
            OutputX = 1. - OutputX;
            _EERepl_OutputX2 = OutputX;
        }
        NofXd2 = OutputX;
        _t49 = strike;
        _t51 = -rate;
        _t50 = time;
        _t52 = _t51 * _t50;
        _t48 = exp(_t52);
        FutureValueX = _t49 * _t48;
        _EERepl_FutureValueX1 = FutureValueX;
        _cond4 = otype == 0;
        if (_cond4)
        {
            _t54 = sptprice;
            _t53 = NofXd1;
            _t56 = FutureValueX;
            _t55 = NofXd2;
            OptionPrice = (_t54 * _t53) - (_t56 * _t55);
            _EERepl_OptionPrice1 = OptionPrice;
        }
        else
        {
            NegNofXd1 = (1. - NofXd1);
            _EERepl_NegNofXd11 = NegNofXd1;
            NegNofXd2 = (1. - NofXd2);
            _EERepl_NegNofXd21 = NegNofXd2;
            _t58 = FutureValueX;
            _t57 = NegNofXd2;
            _t60 = sptprice;
            _t59 = NegNofXd1;
            OptionPrice = (_t58 * _t57) - (_t60 * _t59);
            _EERepl_OptionPrice2 = OptionPrice;
        }
        float BlkSchlsEqEuroNoDiv_return = OptionPrice;
        goto _label0;
    _label0:
        _d_OptionPrice += 1;
        if (_cond4)
        {
            {
                float _r_d77 = _d_OptionPrice;
                float _r70 = _r_d77 * _t53;
                *_d_sptprice += _r70;
                float _r71 = _t54 * _r_d77;
                _d_NofXd1 += _r71;
                float _r72 = -_r_d77 * _t55;
                _d_FutureValueX += _r72;
                float _r73 = _t56 * -_r_d77;
                _d_NofXd2 += _r73;
                _delta_OptionPrice += clad::getErrorVal(_r_d77, _EERepl_OptionPrice1, "OptionPrice");
                _d_OptionPrice -= _r_d77;
            }
        }
        else
        {
            {
                float _r_d80 = _d_OptionPrice;
                float _r74 = _r_d80 * _t57;
                _d_FutureValueX += _r74;
                float _r75 = _t58 * _r_d80;
                _d_NegNofXd2 += _r75;
                float _r76 = -_r_d80 * _t59;
                *_d_sptprice += _r76;
                float _r77 = _t60 * -_r_d80;
                _d_NegNofXd1 += _r77;
                _delta_OptionPrice += clad::getErrorVal(_r_d80, _EERepl_OptionPrice2, "OptionPrice");
                _d_OptionPrice -= _r_d80;
            }
            {
                float _r_d79 = _d_NegNofXd2;
                _d_NofXd2 += -_r_d79;
                _delta_NegNofXd2 += clad::getErrorVal(_r_d79, _EERepl_NegNofXd21, "NegNofXd2");
                _d_NegNofXd2 -= _r_d79;
            }
            {
                float _r_d78 = _d_NegNofXd1;
                _d_NofXd1 += -_r_d78;
                _delta_NegNofXd1 += clad::getErrorVal(_r_d78, _EERepl_NegNofXd11, "NegNofXd1");
                _d_NegNofXd1 -= _r_d78;
            }
        }
        {
            float _r_d76 = _d_FutureValueX;
            double _r65 = _r_d76 * _t48;
            *_d_strike += _r65;
            float _r66 = _t49 * _r_d76;
            double _r67 = _r66 * clad::custom_derivatives::exp_pushforward(_t52, 1.F).pushforward;
            double _r68 = _r67 * _t50;
            *_d_rate += -_r68;
            double _r69 = _t51 * _r67;
            *_d_time += _r69;
            _delta_FutureValueX += clad::getErrorVal(_r_d76, _EERepl_FutureValueX1, "FutureValueX");
            _d_FutureValueX -= _r_d76;
        }
        {
            float _r_d75 = _d_NofXd2;
            _d_OutputX += _r_d75;
            _d_NofXd2 -= _r_d75;
        }
        if (_cond3)
        {
            {
                float _r_d74 = _d_OutputX;
                _d_OutputX += -_r_d74;
                _delta_OutputX += clad::getErrorVal(_r_d74, _EERepl_OutputX2, "OutputX");
                _d_OutputX -= _r_d74;
            }
        }
        {
            float _r_d73 = _d_OutputX;
            _d_xLocal += _r_d73;
            _d_OutputX -= _r_d73;
        }
        {
            float _r_d72 = _d_xLocal;
            _d_xLocal += -_r_d72;
            _delta_xLocal += clad::getErrorVal(_r_d72, _EERepl_xLocal5, "xLocal");
            _d_xLocal -= _r_d72;
        }
        {
            float _r_d71 = _d_xLocal;
            float _r63 = _r_d71 * _t46;
            _d_xLocal_1 += _r63;
            float _r64 = _t47 * _r_d71;
            _d_xNPrimeofX += _r64;
            _delta_xLocal += clad::getErrorVal(_r_d71, _EERepl_xLocal4, "xLocal");
            _d_xLocal -= _r_d71;
        }
        {
            float _r_d70 = _d_xLocal_1;
            _d_xLocal_2 += _r_d70;
            _d_xLocal_1 += _r_d70;
            _delta_xLocal_1 += clad::getErrorVal(_r_d70, _EERepl_xLocal_14, "xLocal_1");
            _d_xLocal_1 -= _r_d70;
        }
        {
            float _r_d69 = _d_xLocal_2;
            _d_xLocal_2 += _r_d69;
            _d_xLocal_3 += _r_d69;
            _delta_xLocal_2 += clad::getErrorVal(_r_d69, _EERepl_xLocal_29, "xLocal_2");
            _d_xLocal_2 -= _r_d69;
        }
        {
            float _r_d68 = _d_xLocal_3;
            double _r62 = _r_d68 * 1.3302744289999999;
            _d_xK2_5 += _r62;
            _delta_xLocal_3 += clad::getErrorVal(_r_d68, _EERepl_xLocal_36, "xLocal_3");
            _d_xLocal_3 -= _r_d68;
        }
        {
            float _r_d67 = _d_xLocal_2;
            _d_xLocal_2 += _r_d67;
            _d_xLocal_3 += _r_d67;
            _delta_xLocal_2 += clad::getErrorVal(_r_d67, _EERepl_xLocal_28, "xLocal_2");
            _d_xLocal_2 -= _r_d67;
        }
        {
            float _r_d66 = _d_xLocal_3;
            double _r61 = _r_d66 * (-1.8212559779999999);
            _d_xK2_4 += _r61;
            _delta_xLocal_3 += clad::getErrorVal(_r_d66, _EERepl_xLocal_35, "xLocal_3");
            _d_xLocal_3 -= _r_d66;
        }
        {
            float _r_d65 = _d_xLocal_2;
            _d_xLocal_2 += _r_d65;
            _d_xLocal_3 += _r_d65;
            _delta_xLocal_2 += clad::getErrorVal(_r_d65, _EERepl_xLocal_27, "xLocal_2");
            _d_xLocal_2 -= _r_d65;
        }
        {
            float _r_d64 = _d_xLocal_3;
            double _r60 = _r_d64 * 1.781477937;
            _d_xK2_3 += _r60;
            _delta_xLocal_3 += clad::getErrorVal(_r_d64, _EERepl_xLocal_34, "xLocal_3");
            _d_xLocal_3 -= _r_d64;
        }
        {
            float _r_d63 = _d_xLocal_2;
            double _r59 = _r_d63 * (-0.356563782);
            _d_xK2_2 += _r59;
            _delta_xLocal_2 += clad::getErrorVal(_r_d63, _EERepl_xLocal_26, "xLocal_2");
            _d_xLocal_2 -= _r_d63;
        }
        {
            float _r_d62 = _d_xLocal_1;
            double _r58 = _r_d62 * 0.31938153000000002;
            _d_xK2 += _r58;
            _delta_xLocal_1 += clad::getErrorVal(_r_d62, _EERepl_xLocal_13, "xLocal_1");
            _d_xLocal_1 -= _r_d62;
        }
        {
            float _r_d61 = _d_xK2_5;
            float _r56 = _r_d61 * _t44;
            _d_xK2_4 += _r56;
            float _r57 = _t45 * _r_d61;
            _d_xK2 += _r57;
            _delta_xK2_5 += clad::getErrorVal(_r_d61, _EERepl_xK2_52, "xK2_5");
            _d_xK2_5 -= _r_d61;
        }
        {
            float _r_d60 = _d_xK2_4;
            float _r54 = _r_d60 * _t42;
            _d_xK2_3 += _r54;
            float _r55 = _t43 * _r_d60;
            _d_xK2 += _r55;
            _delta_xK2_4 += clad::getErrorVal(_r_d60, _EERepl_xK2_43, "xK2_4");
            _d_xK2_4 -= _r_d60;
        }
        {
            float _r_d59 = _d_xK2_3;
            float _r52 = _r_d59 * _t40;
            _d_xK2_2 += _r52;
            float _r53 = _t41 * _r_d59;
            _d_xK2 += _r53;
            _delta_xK2_3 += clad::getErrorVal(_r_d59, _EERepl_xK2_32, "xK2_3");
            _d_xK2_3 -= _r_d59;
        }
        {
            float _r_d58 = _d_xK2_2;
            float _r50 = _r_d58 * _t38;
            _d_xK2 += _r50;
            float _r51 = _t39 * _r_d58;
            _d_xK2 += _r51;
            _delta_xK2_2 += clad::getErrorVal(_r_d58, _EERepl_xK2_23, "xK2_2");
            _d_xK2_2 -= _r_d58;
        }
        {
            float _r_d57 = _d_xK2;
            double _r48 = _r_d57 / _t37;
            double _r49 = _r_d57 * -1. / (_t37 * _t37);
            _d_xK2 += _r49;
            _delta_xK2 += clad::getErrorVal(_r_d57, _EERepl_xK26, "xK2");
            _d_xK2 -= _r_d57;
        }
        {
            float _r_d56 = _d_xK2;
            _d_xK2 += _r_d56;
            _delta_xK2 += clad::getErrorVal(_r_d56, _EERepl_xK25, "xK2");
            _d_xK2 -= _r_d56;
        }
        {
            float _r_d55 = _d_xK2;
            double _r46 = _r_d55 * _t36;
            double _r47 = 0.23164190000000001 * _r_d55;
            _d_xInput += _r47;
            _delta_xK2 += clad::getErrorVal(_r_d55, _EERepl_xK24, "xK2");
            _d_xK2 -= _r_d55;
        }
        {
            float _r_d54 = _d_xNPrimeofX;
            float _r44 = _r_d54 * _t34;
            _d_xNPrimeofX += _r44;
            float _r45 = _t35 * _r_d54;
            _d_inv_sqrt_2xPI += _r45;
            _delta_xNPrimeofX += clad::getErrorVal(_r_d54, _EERepl_xNPrimeofX2, "xNPrimeofX");
            _d_xNPrimeofX -= _r_d54;
        }
        {
            float _r_d53 = _d_xNPrimeofX;
            _d_expValues += _r_d53;
            _d_xNPrimeofX -= _r_d53;
        }
        {
            float _r_d52 = _d_expValues;
            double _r39 = _r_d52 * clad::custom_derivatives::exp_pushforward(_t33, 1.F).pushforward;
            double _r40 = _r39 * _t30;
            double _r41 = _r40 * _t31;
            double _r42 = -0.5F * _r40;
            _d_InputX += _r42;
            double _r43 = _t32 * _r39;
            _d_InputX += _r43;
            _delta_expValues += clad::getErrorVal(_r_d52, _EERepl_expValues2, "expValues");
            _d_expValues -= _r_d52;
        }
        {
            float _r_d51 = _d_xInput;
            _d_InputX += _r_d51;
            _d_xInput -= _r_d51;
        }
        if (_cond2)
        {
            {
                int _r_d49 = _d_sign;
                _d_sign -= _r_d49;
            }
            {
                float _r_d48 = _d_InputX;
                _d_InputX += -_r_d48;
                _delta_InputX += clad::getErrorVal(_r_d48, _EERepl_InputX2, "InputX");
                _d_InputX -= _r_d48;
            }
        }
        else
        {
            int _r_d50 = _d_sign;
            _d_sign -= _r_d50;
        }
        {
            float _r_d47 = _d_InputX;
            _d_d2 += _r_d47;
            _d_InputX -= _r_d47;
        }
        {
            float _r_d46 = _d_NofXd1;
            _d_OutputX += _r_d46;
            _d_NofXd1 -= _r_d46;
        }
        if (_cond1)
        {
            {
                float _r_d45 = _d_OutputX;
                _d_OutputX += -_r_d45;
                _delta_OutputX += clad::getErrorVal(_r_d45, _EERepl_OutputX1, "OutputX");
                _d_OutputX -= _r_d45;
            }
        }
        {
            float _r_d44 = _d_OutputX;
            _d_xLocal += _r_d44;
            _d_OutputX -= _r_d44;
        }
        {
            float _r_d43 = _d_xLocal;
            _d_xLocal += -_r_d43;
            _delta_xLocal += clad::getErrorVal(_r_d43, _EERepl_xLocal3, "xLocal");
            _d_xLocal -= _r_d43;
        }
        {
            float _r_d42 = _d_xLocal;
            float _r37 = _r_d42 * _t28;
            _d_xLocal_1 += _r37;
            float _r38 = _t29 * _r_d42;
            _d_xNPrimeofX += _r38;
            _delta_xLocal += clad::getErrorVal(_r_d42, _EERepl_xLocal2, "xLocal");
            _d_xLocal -= _r_d42;
        }
        {
            float _r_d41 = _d_xLocal_1;
            _d_xLocal_2 += _r_d41;
            _d_xLocal_1 += _r_d41;
            _delta_xLocal_1 += clad::getErrorVal(_r_d41, _EERepl_xLocal_12, "xLocal_1");
            _d_xLocal_1 -= _r_d41;
        }
        {
            float _r_d40 = _d_xLocal_2;
            _d_xLocal_2 += _r_d40;
            _d_xLocal_3 += _r_d40;
            _delta_xLocal_2 += clad::getErrorVal(_r_d40, _EERepl_xLocal_25, "xLocal_2");
            _d_xLocal_2 -= _r_d40;
        }
        {
            float _r_d39 = _d_xLocal_3;
            double _r36 = _r_d39 * 1.3302744289999999;
            _d_xK2_5 += _r36;
            _delta_xLocal_3 += clad::getErrorVal(_r_d39, _EERepl_xLocal_33, "xLocal_3");
            _d_xLocal_3 -= _r_d39;
        }
        {
            float _r_d38 = _d_xLocal_2;
            _d_xLocal_2 += _r_d38;
            _d_xLocal_3 += _r_d38;
            _delta_xLocal_2 += clad::getErrorVal(_r_d38, _EERepl_xLocal_24, "xLocal_2");
            _d_xLocal_2 -= _r_d38;
        }
        {
            float _r_d37 = _d_xLocal_3;
            double _r35 = _r_d37 * (-1.8212559779999999);
            _d_xK2_4 += _r35;
            _delta_xLocal_3 += clad::getErrorVal(_r_d37, _EERepl_xLocal_32, "xLocal_3");
            _d_xLocal_3 -= _r_d37;
        }
        {
            float _r_d36 = _d_xLocal_2;
            _d_xLocal_2 += _r_d36;
            _d_xLocal_3 += _r_d36;
            _delta_xLocal_2 += clad::getErrorVal(_r_d36, _EERepl_xLocal_23, "xLocal_2");
            _d_xLocal_2 -= _r_d36;
        }
        {
            float _r_d35 = _d_xLocal_3;
            double _r34 = _r_d35 * 1.781477937;
            _d_xK2_3 += _r34;
            _delta_xLocal_3 += clad::getErrorVal(_r_d35, _EERepl_xLocal_31, "xLocal_3");
            _d_xLocal_3 -= _r_d35;
        }
        {
            float _r_d34 = _d_xLocal_2;
            double _r33 = _r_d34 * (-0.356563782);
            _d_xK2_2 += _r33;
            _delta_xLocal_2 += clad::getErrorVal(_r_d34, _EERepl_xLocal_22, "xLocal_2");
            _d_xLocal_2 -= _r_d34;
        }
        {
            float _r_d33 = _d_xLocal_1;
            double _r32 = _r_d33 * 0.31938153000000002;
            _d_xK2 += _r32;
            _delta_xLocal_1 += clad::getErrorVal(_r_d33, _EERepl_xLocal_11, "xLocal_1");
            _d_xLocal_1 -= _r_d33;
        }
        {
            float _r_d32 = _d_xK2_5;
            float _r30 = _r_d32 * _t26;
            _d_xK2_4 += _r30;
            float _r31 = _t27 * _r_d32;
            _d_xK2 += _r31;
            _delta_xK2_5 += clad::getErrorVal(_r_d32, _EERepl_xK2_51, "xK2_5");
            _d_xK2_5 -= _r_d32;
        }
        {
            float _r_d31 = _d_xK2_4;
            float _r28 = _r_d31 * _t24;
            _d_xK2_3 += _r28;
            float _r29 = _t25 * _r_d31;
            _d_xK2 += _r29;
            _delta_xK2_4 += clad::getErrorVal(_r_d31, _EERepl_xK2_42, "xK2_4");
            _d_xK2_4 -= _r_d31;
        }
        {
            float _r_d30 = _d_xK2_3;
            float _r26 = _r_d30 * _t22;
            _d_xK2_2 += _r26;
            float _r27 = _t23 * _r_d30;
            _d_xK2 += _r27;
            _delta_xK2_3 += clad::getErrorVal(_r_d30, _EERepl_xK2_31, "xK2_3");
            _d_xK2_3 -= _r_d30;
        }
        {
            float _r_d29 = _d_xK2_2;
            float _r24 = _r_d29 * _t20;
            _d_xK2 += _r24;
            float _r25 = _t21 * _r_d29;
            _d_xK2 += _r25;
            _delta_xK2_2 += clad::getErrorVal(_r_d29, _EERepl_xK2_22, "xK2_2");
            _d_xK2_2 -= _r_d29;
        }
        {
            float _r_d28 = _d_xK2;
            double _r22 = _r_d28 / _t19;
            double _r23 = _r_d28 * -1. / (_t19 * _t19);
            _d_xK2 += _r23;
            _delta_xK2 += clad::getErrorVal(_r_d28, _EERepl_xK23, "xK2");
            _d_xK2 -= _r_d28;
        }
        {
            float _r_d27 = _d_xK2;
            _d_xK2 += _r_d27;
            _delta_xK2 += clad::getErrorVal(_r_d27, _EERepl_xK22, "xK2");
            _d_xK2 -= _r_d27;
        }
        {
            float _r_d26 = _d_xK2;
            double _r20 = _r_d26 * _t18;
            double _r21 = 0.23164190000000001 * _r_d26;
            _d_xInput += _r21;
            _delta_xK2 += clad::getErrorVal(_r_d26, _EERepl_xK21, "xK2");
            _d_xK2 -= _r_d26;
        }
        {
            float _r_d25 = _d_xNPrimeofX;
            float _r18 = _r_d25 * _t16;
            _d_xNPrimeofX += _r18;
            float _r19 = _t17 * _r_d25;
            _d_inv_sqrt_2xPI += _r19;
            _delta_xNPrimeofX += clad::getErrorVal(_r_d25, _EERepl_xNPrimeofX1, "xNPrimeofX");
            _d_xNPrimeofX -= _r_d25;
        }
        {
            float _r_d24 = _d_xNPrimeofX;
            _d_expValues += _r_d24;
            _d_xNPrimeofX -= _r_d24;
        }
        {
            float _r_d23 = _d_expValues;
            double _r13 = _r_d23 * clad::custom_derivatives::exp_pushforward(_t15, 1.F).pushforward;
            double _r14 = _r13 * _t12;
            double _r15 = _r14 * _t13;
            double _r16 = -0.5F * _r14;
            _d_InputX += _r16;
            double _r17 = _t14 * _r13;
            _d_InputX += _r17;
            _delta_expValues += clad::getErrorVal(_r_d23, _EERepl_expValues1, "expValues");
            _d_expValues -= _r_d23;
        }
        {
            float _r_d22 = _d_xInput;
            _d_InputX += _r_d22;
            _d_xInput -= _r_d22;
        }
        if (_cond0)
        {
            {
                int _r_d20 = _d_sign;
                _d_sign -= _r_d20;
            }
            {
                float _r_d19 = _d_InputX;
                _d_InputX += -_r_d19;
                _delta_InputX += clad::getErrorVal(_r_d19, _EERepl_InputX1, "InputX");
                _d_InputX -= _r_d19;
            }
        }
        else
        {
            int _r_d21 = _d_sign;
            _d_sign -= _r_d21;
        }
        {
            float _r_d18 = _d_InputX;
            _d_d1 += _r_d18;
            _d_InputX -= _r_d18;
        }
        {
            float _r_d17 = _d_d2;
            _d_xD2 += _r_d17;
            _d_d2 -= _r_d17;
        }
        {
            float _r_d16 = _d_d1;
            _d_xD1 += _r_d16;
            _d_d1 -= _r_d16;
        }
        {
            float _r_d15 = _d_xD2;
            _d_xD1 += _r_d15;
            _d_xDen += -_r_d15;
            _delta_xD2 += clad::getErrorVal(_r_d15, _EERepl_xD21, "xD2");
            _d_xD2 -= _r_d15;
        }
        {
            float _r_d14 = _d_xD1;
            float _r11 = _r_d14 / _t10;
            _d_xD1 += _r11;
            float _r12 = _r_d14 * -_t11 / (_t10 * _t10);
            _d_xDen += _r12;
            _delta_xD1 += clad::getErrorVal(_r_d14, _EERepl_xD14, "xD1");
            _d_xD1 -= _r_d14;
        }
        {
            float _r_d13 = _d_xDen;
            float _r9 = _r_d13 * _t8;
            _d_xVolatility += _r9;
            float _r10 = _t9 * _r_d13;
            _d_xSqrtTime += _r10;
            _delta_xDen += clad::getErrorVal(_r_d13, _EERepl_xDen1, "xDen");
            _d_xDen -= _r_d13;
        }
        {
            float _r_d12 = _d_xD1;
            _d_xD1 += _r_d12;
            _d_xLogTerm += _r_d12;
            _delta_xD1 += clad::getErrorVal(_r_d12, _EERepl_xD13, "xD1");
            _d_xD1 -= _r_d12;
        }
        {
            float _r_d11 = _d_xD1;
            float _r7 = _r_d11 * _t6;
            _d_xD1 += _r7;
            float _r8 = _t7 * _r_d11;
            _d_xTime += _r8;
            _delta_xD1 += clad::getErrorVal(_r_d11, _EERepl_xD12, "xD1");
            _d_xD1 -= _r_d11;
        }
        {
            float _r_d10 = _d_xD1;
            _d_xRiskFreeRate += _r_d10;
            _d_xPowerTerm += _r_d10;
            _delta_xD1 += clad::getErrorVal(_r_d10, _EERepl_xD11, "xD1");
            _d_xD1 -= _r_d10;
        }
        {
            float _r_d9 = _d_xPowerTerm;
            double _r6 = _r_d9 * 0.5;
            _d_xPowerTerm += _r6;
            _delta_xPowerTerm += clad::getErrorVal(_r_d9, _EERepl_xPowerTerm2, "xPowerTerm");
            _d_xPowerTerm -= _r_d9;
        }
        {
            float _r_d8 = _d_xPowerTerm;
            float _r4 = _r_d8 * _t4;
            _d_xVolatility += _r4;
            float _r5 = _t5 * _r_d8;
            _d_xVolatility += _r5;
            _delta_xPowerTerm += clad::getErrorVal(_r_d8, _EERepl_xPowerTerm1, "xPowerTerm");
            _d_xPowerTerm -= _r_d8;
        }
        {
            float _r_d7 = _d_xLogTerm;
            _d_logValues += _r_d7;
            _d_xLogTerm -= _r_d7;
        }
        {
            float _r_d6 = _d_logValues;
            double _r1 = _r_d6 * clad::custom_derivatives::log_pushforward(_t3, 1.F).pushforward;
            double _r2 = _r1 / _t1;
            *_d_sptprice += _r2;
            double _r3 = _r1 * -_t2 / (_t1 * _t1);
            *_d_strike += _r3;
            _delta_logValues += clad::getErrorVal(_r_d6, _EERepl_logValues1, "logValues");
            _d_logValues -= _r_d6;
        }
        {
            float _r_d5 = _d_xSqrtTime;
            double _r0 = _r_d5 * clad::custom_derivatives::sqrt_pushforward(_t0, 1.F).pushforward;
            _d_xTime += _r0;
            _delta_xSqrtTime += clad::getErrorVal(_r_d5, _EERepl_xSqrtTime1, "xSqrtTime");
            _d_xSqrtTime -= _r_d5;
        }
        {
            float _r_d4 = _d_xTime;
            *_d_time += _r_d4;
            _d_xTime -= _r_d4;
        }
        {
            float _r_d3 = _d_xVolatility;
            *_d_volatility += _r_d3;
            _d_xVolatility -= _r_d3;
        }
        {
            float _r_d2 = _d_xRiskFreeRate;
            *_d_rate += _r_d2;
            _d_xRiskFreeRate -= _r_d2;
        }
        {
            float _r_d1 = _d_xStrikePrice;
            *_d_strike += _r_d1;
            _d_xStrikePrice -= _r_d1;
        }
        {
            float _r_d0 = _d_xStockPrice;
            *_d_sptprice += _r_d0;
            _d_xStockPrice -= _r_d0;
        }
        _delta_inv_sqrt_2xPI += clad::getErrorVal(_d_inv_sqrt_2xPI, _EERepl_inv_sqrt_2xPI0, "inv_sqrt_2xPI");
        double _delta_sptprice = 0;
        _delta_sptprice += clad::getErrorVal(*_d_sptprice, sptprice, "sptprice");
        double _delta_strike = 0;
        _delta_strike += clad::getErrorVal(*_d_strike, strike, "strike");
        double _delta_rate = 0;
        _delta_rate += clad::getErrorVal(*_d_rate, rate, "rate");
        double _delta_volatility = 0;
        _delta_volatility += clad::getErrorVal(*_d_volatility, volatility, "volatility");
        double _delta_time = 0;
        _delta_time += clad::getErrorVal(*_d_time, time, "time");
        _final_error += _delta_time + _delta_volatility + _delta_xLocal_3 + _delta_xK2_5 + _delta_xK2_3 + _delta_xLocal_2 + _delta_xDen + _delta_xPowerTerm + _delta_sptprice + _delta_inv_sqrt_2xPI + _delta_xD2 + _delta_xD1 + _delta_xLogTerm + _delta_InputX + _delta_NegNofXd2 + _delta_xStockPrice + _delta_NegNofXd1 + _delta_rate + _delta_OptionPrice + _delta_NofXd2 + _delta_xLocal_1 + _delta_d1 + _delta_d2 + _delta_FutureValueX + _delta_NofXd1 + _delta_xStrikePrice + _delta_OutputX + _delta_xRiskFreeRate + _delta_xInput + _delta_xVolatility + _delta_xNPrimeofX + _delta_strike + _delta_xTime + _delta_expValues + _delta_xSqrtTime + _delta_xLocal + _delta_xK2 + _delta_xK2_4 + _delta_logValues + _delta_xK2_2;
    }

}
