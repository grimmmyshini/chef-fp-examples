
namespace approx {
void BlkSchlsEqEuroNoDiv_grad(double sptprice, double strike, double rate, double volatility, double time, int otype, clad::array_ref<double> _d_sptprice, clad::array_ref<double> _d_strike, clad::array_ref<double> _d_rate, clad::array_ref<double> _d_volatility, clad::array_ref<double> _d_time, clad::array_ref<int> _d_otype, double &_final_error) {
    double _delta_OptionPrice = 0;
    double _EERepl_OptionPrice0;
    double _delta_xStockPrice = 0;
    double _EERepl_xStockPrice0;
    double _delta_xStrikePrice = 0;
    double _EERepl_xStrikePrice0;
    double _delta_xRiskFreeRate = 0;
    double _EERepl_xRiskFreeRate0;
    double _delta_xVolatility = 0;
    double _EERepl_xVolatility0;
    double _delta_xTime = 0;
    double _EERepl_xTime0;
    double _delta_xSqrtTime = 0;
    double _EERepl_xSqrtTime0;
    double _delta_logValues = 0;
    double _EERepl_logValues0;
    double _delta_xLogTerm = 0;
    double _EERepl_xLogTerm0;
    double _delta_xD1 = 0;
    double _EERepl_xD10;
    double _delta_xD2 = 0;
    double _EERepl_xD20;
    double _delta_xPowerTerm = 0;
    double _EERepl_xPowerTerm0;
    double _delta_xDen = 0;
    double _EERepl_xDen0;
    double _delta_d1 = 0;
    double _EERepl_d10;
    double _delta_d2 = 0;
    double _EERepl_d20;
    double _delta_FutureValueX = 0;
    double _EERepl_FutureValueX0;
    double _delta_NofXd1 = 0;
    double _EERepl_NofXd10;
    double _delta_NofXd2 = 0;
    double _EERepl_NofXd20;
    double _delta_NegNofXd1 = 0;
    double _EERepl_NegNofXd10;
    double _delta_NegNofXd2 = 0;
    double _EERepl_NegNofXd20;
    double _delta_InputX = 0;
    double _EERepl_InputX0;
    double _delta_OutputX = 0;
    double _EERepl_OutputX0;
    double _delta_xInput = 0;
    double _EERepl_xInput0;
    double _delta_xNPrimeofX = 0;
    double _EERepl_xNPrimeofX0;
    double _delta_expValues = 0;
    double _EERepl_expValues0;
    double _delta_xK2 = 0;
    double _EERepl_xK20;
    double _delta_xK2_2 = 0;
    double _EERepl_xK2_20;
    double _delta_xK2_20 = 0;
    double _EERepl_xK2_21;
    double _delta_xK2_4 = 0;
    double _EERepl_xK2_40;
    double _delta_xK2_40 = 0;
    double _EERepl_xK2_41;
    double _delta_xLocal = 0;
    double _EERepl_xLocal0;
    double _delta_xLocal0 = 0;
    double _EERepl_xLocal1;
    double _delta_xLocal_2 = 0;
    double _EERepl_xLocal_20;
    double _delta_xLocal_20 = 0;
    double _EERepl_xLocal_21;
    double _delta_inv_sqrt_2xPI = 0;
    double _EERepl_inv_sqrt_2xPI0;
    double _delta_clad_log_logValues_ = 0;
    double _EERepl_clad_log_logValues_0;
    double _delta_clad_sqr_xSqrtTime_ = 0;
    double _EERepl_clad_sqr_xSqrtTime_0;
    double _delta_clad_exp_expValues_ = 0;
    double _EERepl_clad_exp_expValues_0;
    double _delta_clad_exp_FutureValueX_ = 0;
    double _EERepl_clad_exp_FutureValueX_0;
    double _t0;
    double _EERepl_xSqrtTime1;
    double _t1;
    double _t2;
    double _EERepl_clad_log_logValues_1;
    double _t3;
    double _EERepl_logValues1;
    double _t4;
    double _t5;
    double _EERepl_xPowerTerm1;
    double _EERepl_xPowerTerm2;
    double _EERepl_xD11;
    double _t6;
    double _t7;
    double _EERepl_xD12;
    double _EERepl_xD13;
    double _t8;
    double _t9;
    double _EERepl_xDen1;
    double _t10;
    double _t11;
    double _EERepl_xD14;
    double _EERepl_xD21;
    bool _cond0;
    double _EERepl_InputX1;
    double _t12;
    double _t13;
    double _t14;
    double _EERepl_clad_exp_expValues_1;
    double _t15;
    double _EERepl_expValues1;
    double _t16;
    double _t17;
    double _EERepl_xNPrimeofX1;
    double _t18;
    double _EERepl_xK21;
    double _EERepl_xK22;
    double _t19;
    double _EERepl_xK23;
    double _t20;
    double _t21;
    double _EERepl_xK2_22;
    double _delta_xK2_3 = 0;
    double _t22;
    double _t23;
    double _EERepl_xK2_31;
    double _t24;
    double _t25;
    double _EERepl_xK2_42;
    double _delta_xK2_5 = 0;
    double _t26;
    double _t27;
    double _EERepl_xK2_51;
    double _delta_xLocal_1 = 0;
    double _EERepl_xLocal_11;
    double _EERepl_xLocal_22;
    double _delta_xLocal_3 = 0;
    double _EERepl_xLocal_31;
    double _EERepl_xLocal_23;
    double _EERepl_xLocal_32;
    double _EERepl_xLocal_24;
    double _EERepl_xLocal_33;
    double _EERepl_xLocal_25;
    double _EERepl_xLocal_12;
    double _t28;
    double _t29;
    double _EERepl_xLocal2;
    double _EERepl_xLocal3;
    int _cond1;
    double _EERepl_OutputX1;
    bool _cond2;
    double _EERepl_InputX2;
    double _t30;
    double _t31;
    double _t32;
    double _EERepl_clad_exp_expValues_2;
    double _t33;
    double _EERepl_expValues2;
    double _t34;
    double _t35;
    double _EERepl_xNPrimeofX2;
    double _t36;
    double _EERepl_xK24;
    double _EERepl_xK25;
    double _t37;
    double _EERepl_xK26;
    double _t38;
    double _t39;
    double _EERepl_xK2_23;
    double _t40;
    double _t41;
    double _EERepl_xK2_32;
    double _t42;
    double _t43;
    double _EERepl_xK2_43;
    double _t44;
    double _t45;
    double _EERepl_xK2_52;
    double _EERepl_xLocal_13;
    double _EERepl_xLocal_26;
    double _EERepl_xLocal_34;
    double _EERepl_xLocal_27;
    double _EERepl_xLocal_35;
    double _EERepl_xLocal_28;
    double _EERepl_xLocal_36;
    double _EERepl_xLocal_29;
    double _EERepl_xLocal_14;
    double _t46;
    double _t47;
    double _EERepl_xLocal4;
    double _EERepl_xLocal5;
    int _cond3;
    double _EERepl_OutputX2;
    double _t48;
    double _t49;
    double _EERepl_clad_exp_FutureValueX_1;
    double _t50;
    double _t51;
    double _t52;
    double _EERepl_FutureValueX1;
    bool _cond4;
    double _t53;
    double _t54;
    double _t55;
    double _t56;
    double _EERepl_OptionPrice1;
    double _EERepl_NegNofXd11;
    double _EERepl_NegNofXd21;
    double _t57;
    double _t58;
    double _t59;
    double _t60;
    double _EERepl_OptionPrice2;
    double _d_OptionPrice = 0;
    double OptionPrice;
    _EERepl_OptionPrice0 = OptionPrice;
    double _d_xStockPrice = 0;
    double xStockPrice;
    _EERepl_xStockPrice0 = xStockPrice;
    double _d_xStrikePrice = 0;
    double xStrikePrice;
    _EERepl_xStrikePrice0 = xStrikePrice;
    double _d_xRiskFreeRate = 0;
    double xRiskFreeRate;
    _EERepl_xRiskFreeRate0 = xRiskFreeRate;
    double _d_xVolatility = 0;
    double xVolatility;
    _EERepl_xVolatility0 = xVolatility;
    double _d_xTime = 0;
    double xTime;
    _EERepl_xTime0 = xTime;
    double _d_xSqrtTime = 0;
    double xSqrtTime;
    _EERepl_xSqrtTime0 = xSqrtTime;
    double _d_logValues = 0;
    double logValues;
    _EERepl_logValues0 = logValues;
    double _d_xLogTerm = 0;
    double xLogTerm;
    _EERepl_xLogTerm0 = xLogTerm;
    double _d_xD1 = 0;
    double xD1;
    _EERepl_xD10 = xD1;
    double _d_xD2 = 0;
    double xD2;
    _EERepl_xD20 = xD2;
    double _d_xPowerTerm = 0;
    double xPowerTerm;
    _EERepl_xPowerTerm0 = xPowerTerm;
    double _d_xDen = 0;
    double xDen;
    _EERepl_xDen0 = xDen;
    double _d_d1 = 0;
    double d1;
    _EERepl_d10 = d1;
    double _d_d2 = 0;
    double d2;
    _EERepl_d20 = d2;
    double _d_FutureValueX = 0;
    double FutureValueX;
    _EERepl_FutureValueX0 = FutureValueX;
    double _d_NofXd1 = 0;
    double NofXd1;
    _EERepl_NofXd10 = NofXd1;
    double _d_NofXd2 = 0;
    double NofXd2;
    _EERepl_NofXd20 = NofXd2;
    double _d_NegNofXd1 = 0;
    double NegNofXd1;
    _EERepl_NegNofXd10 = NegNofXd1;
    double _d_NegNofXd2 = 0;
    double NegNofXd2;
    _EERepl_NegNofXd20 = NegNofXd2;
    double _d_InputX = 0;
    double InputX;
    _EERepl_InputX0 = InputX;
    int _d_sign = 0;
    int sign;
    double _d_OutputX = 0;
    double OutputX;
    _EERepl_OutputX0 = OutputX;
    double _d_xInput = 0;
    double xInput;
    _EERepl_xInput0 = xInput;
    double _d_xNPrimeofX = 0;
    double xNPrimeofX;
    _EERepl_xNPrimeofX0 = xNPrimeofX;
    double _d_expValues = 0;
    double expValues;
    _EERepl_expValues0 = expValues;
    double _d_xK2 = 0;
    double xK2;
    _EERepl_xK20 = xK2;
    double _d_xK2_2 = 0, _d_xK2_3 = 0;
    double xK2_2, xK2_3;
    _EERepl_xK2_21 = xK2_2;
    _EERepl_xK2_20 = xK2_2;
    double _d_xK2_4 = 0, _d_xK2_5 = 0;
    double xK2_4, xK2_5;
    _EERepl_xK2_41 = xK2_4;
    _EERepl_xK2_40 = xK2_4;
    double _d_xLocal = 0, _d_xLocal_1 = 0;
    double xLocal, xLocal_1;
    _EERepl_xLocal1 = xLocal;
    _EERepl_xLocal0 = xLocal;
    double _d_xLocal_2 = 0, _d_xLocal_3 = 0;
    double xLocal_2, xLocal_3;
    _EERepl_xLocal_21 = xLocal_2;
    _EERepl_xLocal_20 = xLocal_2;
    double _d_inv_sqrt_2xPI = 0;
    double inv_sqrt_2xPI = 0.3989422804014327;
    _EERepl_inv_sqrt_2xPI0 = inv_sqrt_2xPI;
    double _d_clad_log_logValues_ = 0;
    double clad_log_logValues_;
    _EERepl_clad_log_logValues_0 = clad_log_logValues_;
    double _d_clad_sqr_xSqrtTime_ = 0;
    double clad_sqr_xSqrtTime_;
    _EERepl_clad_sqr_xSqrtTime_0 = clad_sqr_xSqrtTime_;
    double _d_clad_exp_expValues_ = 0;
    double clad_exp_expValues_;
    _EERepl_clad_exp_expValues_0 = clad_exp_expValues_;
    double _d_clad_exp_FutureValueX_ = 0;
    double clad_exp_FutureValueX_;
    _EERepl_clad_exp_FutureValueX_0 = clad_exp_FutureValueX_;
    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;
    xTime = time;
    clad_sqr_xSqrtTime_ = xTime;
    double _EERepl_clad_sqr_xSqrtTime_ = clad_sqr_xSqrtTime_;
    _t0 = clad_sqr_xSqrtTime_;
    xSqrtTime = sqrt(_t0);
    _EERepl_xSqrtTime1 = xSqrtTime;
    _t2 = sptprice;
    _t1 = strike;
    clad_log_logValues_ = _t2 / _t1;
    _EERepl_clad_log_logValues_1 = clad_log_logValues_;
    _t3 = clad_log_logValues_;
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
    if (_cond0) {
        InputX = -InputX;
        _EERepl_InputX1 = InputX;
        sign = 1;
    } else
        sign = 0;
    xInput = InputX;
    _t13 = InputX;
    _t14 = -0.5F * _t13;
    _t12 = InputX;
    clad_exp_expValues_ = _t14 * _t12;
    _EERepl_clad_exp_expValues_1 = clad_exp_expValues_;
    _t15 = clad_exp_expValues_;
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
    if (_cond1) {
        OutputX = 1. - OutputX;
        _EERepl_OutputX1 = OutputX;
    }
    NofXd1 = OutputX;
    InputX = d2;
    _cond2 = InputX < 0.;
    if (_cond2) {
        InputX = -InputX;
        _EERepl_InputX2 = InputX;
        sign = 1;
    } else
        sign = 0;
    xInput = InputX;
    _t31 = InputX;
    _t32 = -0.5F * _t31;
    _t30 = InputX;
    clad_exp_expValues_ = _t32 * _t30;
    _EERepl_clad_exp_expValues_2 = clad_exp_expValues_;
    _t33 = clad_exp_expValues_;
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
    if (_cond3) {
        OutputX = 1. - OutputX;
        _EERepl_OutputX2 = OutputX;
    }
    NofXd2 = OutputX;
    _t49 = -rate;
    _t48 = time;
    clad_exp_FutureValueX_ = _t49 * _t48;
    _EERepl_clad_exp_FutureValueX_1 = clad_exp_FutureValueX_;
    _t51 = strike;
    _t52 = clad_exp_FutureValueX_;
    _t50 = exp(_t52);
    FutureValueX = _t51 * _t50;
    _EERepl_FutureValueX1 = FutureValueX;
    _cond4 = otype == 0;
    if (_cond4) {
        _t54 = sptprice;
        _t53 = NofXd1;
        _t56 = FutureValueX;
        _t55 = NofXd2;
        OptionPrice = (_t54 * _t53) - (_t56 * _t55);
        _EERepl_OptionPrice1 = OptionPrice;
    } else {
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
    double BlkSchlsEqEuroNoDiv_return = OptionPrice;
    goto _label0;
  _label0:
    _d_OptionPrice += 1;
    if (_cond4) {
        {
            double _r_d82 = _d_OptionPrice;
            double _r70 = _r_d82 * _t53;
            * _d_sptprice += _r70;
            double _r71 = _t54 * _r_d82;
            _d_NofXd1 += _r71;
            double _r72 = -_r_d82 * _t55;
            _d_FutureValueX += _r72;
            double _r73 = _t56 * -_r_d82;
            _d_NofXd2 += _r73;
            _d_OptionPrice -= _r_d82;
        }
    } else {
        {
            double _r_d85 = _d_OptionPrice;
            double _r74 = _r_d85 * _t57;
            _d_FutureValueX += _r74;
            double _r75 = _t58 * _r_d85;
            _d_NegNofXd2 += _r75;
            double _r76 = -_r_d85 * _t59;
            * _d_sptprice += _r76;
            double _r77 = _t60 * -_r_d85;
            _d_NegNofXd1 += _r77;
            _d_OptionPrice -= _r_d85;
        }
        {
            double _r_d84 = _d_NegNofXd2;
            _d_NofXd2 += -_r_d84;
            _d_NegNofXd2 -= _r_d84;
        }
        {
            double _r_d83 = _d_NegNofXd1;
            _d_NofXd1 += -_r_d83;
            _d_NegNofXd1 -= _r_d83;
        }
    }
    {
        double _r_d81 = _d_FutureValueX;
        double _r67 = _r_d81 * _t50;
        * _d_strike += _r67;
        double _r68 = _t51 * _r_d81;
        double _r69 = _r68 * clad::custom_derivatives::exp_pushforward(_t52, 1.).pushforward;
        _d_clad_exp_FutureValueX_ += _r69;
        _d_FutureValueX -= _r_d81;
    }
    {
        double _r_d80 = _d_clad_exp_FutureValueX_;
        double _r65 = _r_d80 * _t48;
        * _d_rate += -_r65;
        double _r66 = _t49 * _r_d80;
        * _d_time += _r66;
        _delta_clad_exp_FutureValueX_ += clad::doApprox(_r_d80, _EERepl_clad_exp_FutureValueX_1, "clad_exp_FutureValueX");
        _d_clad_exp_FutureValueX_ -= _r_d80;
    }
    {
        double _r_d79 = _d_NofXd2;
        _d_OutputX += _r_d79;
        _d_NofXd2 -= _r_d79;
    }
    if (_cond3) {
        {
            double _r_d78 = _d_OutputX;
            _d_OutputX += -_r_d78;
            _d_OutputX -= _r_d78;
        }
    }
    {
        double _r_d77 = _d_OutputX;
        _d_xLocal += _r_d77;
        _d_OutputX -= _r_d77;
    }
    {
        double _r_d76 = _d_xLocal;
        _d_xLocal += -_r_d76;
        _d_xLocal -= _r_d76;
    }
    {
        double _r_d75 = _d_xLocal;
        double _r63 = _r_d75 * _t46;
        _d_xLocal_1 += _r63;
        double _r64 = _t47 * _r_d75;
        _d_xNPrimeofX += _r64;
        _d_xLocal -= _r_d75;
    }
    {
        double _r_d74 = _d_xLocal_1;
        _d_xLocal_2 += _r_d74;
        _d_xLocal_1 += _r_d74;
        _d_xLocal_1 -= _r_d74;
    }
    {
        double _r_d73 = _d_xLocal_2;
        _d_xLocal_2 += _r_d73;
        _d_xLocal_3 += _r_d73;
        _d_xLocal_2 -= _r_d73;
    }
    {
        double _r_d72 = _d_xLocal_3;
        double _r62 = _r_d72 * 1.3302744289999999;
        _d_xK2_5 += _r62;
        _d_xLocal_3 -= _r_d72;
    }
    {
        double _r_d71 = _d_xLocal_2;
        _d_xLocal_2 += _r_d71;
        _d_xLocal_3 += _r_d71;
        _d_xLocal_2 -= _r_d71;
    }
    {
        double _r_d70 = _d_xLocal_3;
        double _r61 = _r_d70 * (-1.8212559779999999);
        _d_xK2_4 += _r61;
        _d_xLocal_3 -= _r_d70;
    }
    {
        double _r_d69 = _d_xLocal_2;
        _d_xLocal_2 += _r_d69;
        _d_xLocal_3 += _r_d69;
        _d_xLocal_2 -= _r_d69;
    }
    {
        double _r_d68 = _d_xLocal_3;
        double _r60 = _r_d68 * 1.781477937;
        _d_xK2_3 += _r60;
        _d_xLocal_3 -= _r_d68;
    }
    {
        double _r_d67 = _d_xLocal_2;
        double _r59 = _r_d67 * (-0.356563782);
        _d_xK2_2 += _r59;
        _d_xLocal_2 -= _r_d67;
    }
    {
        double _r_d66 = _d_xLocal_1;
        double _r58 = _r_d66 * 0.31938153000000002;
        _d_xK2 += _r58;
        _d_xLocal_1 -= _r_d66;
    }
    {
        double _r_d65 = _d_xK2_5;
        double _r56 = _r_d65 * _t44;
        _d_xK2_4 += _r56;
        double _r57 = _t45 * _r_d65;
        _d_xK2 += _r57;
        _d_xK2_5 -= _r_d65;
    }
    {
        double _r_d64 = _d_xK2_4;
        double _r54 = _r_d64 * _t42;
        _d_xK2_3 += _r54;
        double _r55 = _t43 * _r_d64;
        _d_xK2 += _r55;
        _d_xK2_4 -= _r_d64;
    }
    {
        double _r_d63 = _d_xK2_3;
        double _r52 = _r_d63 * _t40;
        _d_xK2_2 += _r52;
        double _r53 = _t41 * _r_d63;
        _d_xK2 += _r53;
        _d_xK2_3 -= _r_d63;
    }
    {
        double _r_d62 = _d_xK2_2;
        double _r50 = _r_d62 * _t38;
        _d_xK2 += _r50;
        double _r51 = _t39 * _r_d62;
        _d_xK2 += _r51;
        _d_xK2_2 -= _r_d62;
    }
    {
        double _r_d61 = _d_xK2;
        double _r48 = _r_d61 / _t37;
        double _r49 = _r_d61 * -1. / (_t37 * _t37);
        _d_xK2 += _r49;
        _d_xK2 -= _r_d61;
    }
    {
        double _r_d60 = _d_xK2;
        _d_xK2 += _r_d60;
        _d_xK2 -= _r_d60;
    }
    {
        double _r_d59 = _d_xK2;
        double _r46 = _r_d59 * _t36;
        double _r47 = 0.23164190000000001 * _r_d59;
        _d_xInput += _r47;
        _d_xK2 -= _r_d59;
    }
    {
        double _r_d58 = _d_xNPrimeofX;
        double _r44 = _r_d58 * _t34;
        _d_xNPrimeofX += _r44;
        double _r45 = _t35 * _r_d58;
        _d_inv_sqrt_2xPI += _r45;
        _d_xNPrimeofX -= _r_d58;
    }
    {
        double _r_d57 = _d_xNPrimeofX;
        _d_expValues += _r_d57;
        _d_xNPrimeofX -= _r_d57;
    }
    {
        double _r_d56 = _d_expValues;
        double _r43 = _r_d56 * clad::custom_derivatives::exp_pushforward(_t33, 1.).pushforward;
        _d_clad_exp_expValues_ += _r43;
        _d_expValues -= _r_d56;
    }
    {
        double _r_d55 = _d_clad_exp_expValues_;
        double _r39 = _r_d55 * _t30;
        double _r40 = _r39 * _t31;
        double _r41 = -0.5F * _r39;
        _d_InputX += _r41;
        double _r42 = _t32 * _r_d55;
        _d_InputX += _r42;
        _delta_clad_exp_expValues_ += clad::doApprox(_r_d55, _EERepl_clad_exp_expValues_2, "clad_exp_expValues_");
        _d_clad_exp_expValues_ -= _r_d55;
    }
    {
        double _r_d54 = _d_xInput;
        _d_InputX += _r_d54;
        _d_xInput -= _r_d54;
    }
    if (_cond2) {
        {
            int _r_d52 = _d_sign;
            _d_sign -= _r_d52;
        }
        {
            double _r_d51 = _d_InputX;
            _d_InputX += -_r_d51;
            _d_InputX -= _r_d51;
        }
    } else {
        int _r_d53 = _d_sign;
        _d_sign -= _r_d53;
    }
    {
        double _r_d50 = _d_InputX;
        _d_d2 += _r_d50;
        _d_InputX -= _r_d50;
    }
    {
        double _r_d49 = _d_NofXd1;
        _d_OutputX += _r_d49;
        _d_NofXd1 -= _r_d49;
    }
    if (_cond1) {
        {
            double _r_d48 = _d_OutputX;
            _d_OutputX += -_r_d48;
            _d_OutputX -= _r_d48;
        }
    }
    {
        double _r_d47 = _d_OutputX;
        _d_xLocal += _r_d47;
        _d_OutputX -= _r_d47;
    }
    {
        double _r_d46 = _d_xLocal;
        _d_xLocal += -_r_d46;
        _d_xLocal -= _r_d46;
    }
    {
        double _r_d45 = _d_xLocal;
        double _r37 = _r_d45 * _t28;
        _d_xLocal_1 += _r37;
        double _r38 = _t29 * _r_d45;
        _d_xNPrimeofX += _r38;
        _d_xLocal -= _r_d45;
    }
    {
        double _r_d44 = _d_xLocal_1;
        _d_xLocal_2 += _r_d44;
        _d_xLocal_1 += _r_d44;
        _d_xLocal_1 -= _r_d44;
    }
    {
        double _r_d43 = _d_xLocal_2;
        _d_xLocal_2 += _r_d43;
        _d_xLocal_3 += _r_d43;
        _d_xLocal_2 -= _r_d43;
    }
    {
        double _r_d42 = _d_xLocal_3;
        double _r36 = _r_d42 * 1.3302744289999999;
        _d_xK2_5 += _r36;
        _d_xLocal_3 -= _r_d42;
    }
    {
        double _r_d41 = _d_xLocal_2;
        _d_xLocal_2 += _r_d41;
        _d_xLocal_3 += _r_d41;
        _d_xLocal_2 -= _r_d41;
    }
    {
        double _r_d40 = _d_xLocal_3;
        double _r35 = _r_d40 * (-1.8212559779999999);
        _d_xK2_4 += _r35;
        _d_xLocal_3 -= _r_d40;
    }
    {
        double _r_d39 = _d_xLocal_2;
        _d_xLocal_2 += _r_d39;
        _d_xLocal_3 += _r_d39;
        _d_xLocal_2 -= _r_d39;
    }
    {
        double _r_d38 = _d_xLocal_3;
        double _r34 = _r_d38 * 1.781477937;
        _d_xK2_3 += _r34;
        _d_xLocal_3 -= _r_d38;
    }
    {
        double _r_d37 = _d_xLocal_2;
        double _r33 = _r_d37 * (-0.356563782);
        _d_xK2_2 += _r33;
        _d_xLocal_2 -= _r_d37;
    }
    {
        double _r_d36 = _d_xLocal_1;
        double _r32 = _r_d36 * 0.31938153000000002;
        _d_xK2 += _r32;
        _d_xLocal_1 -= _r_d36;
    }
    {
        double _r_d35 = _d_xK2_5;
        double _r30 = _r_d35 * _t26;
        _d_xK2_4 += _r30;
        double _r31 = _t27 * _r_d35;
        _d_xK2 += _r31;
        _d_xK2_5 -= _r_d35;
    }
    {
        double _r_d34 = _d_xK2_4;
        double _r28 = _r_d34 * _t24;
        _d_xK2_3 += _r28;
        double _r29 = _t25 * _r_d34;
        _d_xK2 += _r29;
        _d_xK2_4 -= _r_d34;
    }
    {
        double _r_d33 = _d_xK2_3;
        double _r26 = _r_d33 * _t22;
        _d_xK2_2 += _r26;
        double _r27 = _t23 * _r_d33;
        _d_xK2 += _r27;
        _d_xK2_3 -= _r_d33;
    }
    {
        double _r_d32 = _d_xK2_2;
        double _r24 = _r_d32 * _t20;
        _d_xK2 += _r24;
        double _r25 = _t21 * _r_d32;
        _d_xK2 += _r25;
        _d_xK2_2 -= _r_d32;
    }
    {
        double _r_d31 = _d_xK2;
        double _r22 = _r_d31 / _t19;
        double _r23 = _r_d31 * -1. / (_t19 * _t19);
        _d_xK2 += _r23;
        _d_xK2 -= _r_d31;
    }
    {
        double _r_d30 = _d_xK2;
        _d_xK2 += _r_d30;
        _d_xK2 -= _r_d30;
    }
    {
        double _r_d29 = _d_xK2;
        double _r20 = _r_d29 * _t18;
        double _r21 = 0.23164190000000001 * _r_d29;
        _d_xInput += _r21;
        _d_xK2 -= _r_d29;
    }
    {
        double _r_d28 = _d_xNPrimeofX;
        double _r18 = _r_d28 * _t16;
        _d_xNPrimeofX += _r18;
        double _r19 = _t17 * _r_d28;
        _d_inv_sqrt_2xPI += _r19;
        _d_xNPrimeofX -= _r_d28;
    }
    {
        double _r_d27 = _d_xNPrimeofX;
        _d_expValues += _r_d27;
        _d_xNPrimeofX -= _r_d27;
    }
    {
        double _r_d26 = _d_expValues;
        double _r17 = _r_d26 * clad::custom_derivatives::exp_pushforward(_t15, 1.).pushforward;
        _d_clad_exp_expValues_ += _r17;
        _d_expValues -= _r_d26;
    }
    {
        double _r_d25 = _d_clad_exp_expValues_;
        double _r13 = _r_d25 * _t12;
        double _r14 = _r13 * _t13;
        double _r15 = -0.5F * _r13;
        _d_InputX += _r15;
        double _r16 = _t14 * _r_d25;
        _d_InputX += _r16;
        _delta_clad_exp_expValues_ += clad::doApprox(_r_d25, _EERepl_clad_exp_expValues_1, "clad_exp_expValues_");
        _d_clad_exp_expValues_ -= _r_d25;
    }
    {
        double _r_d24 = _d_xInput;
        _d_InputX += _r_d24;
        _d_xInput -= _r_d24;
    }
    if (_cond0) {
        {
            int _r_d22 = _d_sign;
            _d_sign -= _r_d22;
        }
        {
            double _r_d21 = _d_InputX;
            _d_InputX += -_r_d21;
            _d_InputX -= _r_d21;
        }
    } else {
        int _r_d23 = _d_sign;
        _d_sign -= _r_d23;
    }
    {
        double _r_d20 = _d_InputX;
        _d_d1 += _r_d20;
        _d_InputX -= _r_d20;
    }
    {
        double _r_d19 = _d_d2;
        _d_xD2 += _r_d19;
        _d_d2 -= _r_d19;
    }
    {
        double _r_d18 = _d_d1;
        _d_xD1 += _r_d18;
        _d_d1 -= _r_d18;
    }
    {
        double _r_d17 = _d_xD2;
        _d_xD1 += _r_d17;
        _d_xDen += -_r_d17;
        _d_xD2 -= _r_d17;
    }
    {
        double _r_d16 = _d_xD1;
        double _r11 = _r_d16 / _t10;
        _d_xD1 += _r11;
        double _r12 = _r_d16 * -_t11 / (_t10 * _t10);
        _d_xDen += _r12;
        _d_xD1 -= _r_d16;
    }
    {
        double _r_d15 = _d_xDen;
        double _r9 = _r_d15 * _t8;
        _d_xVolatility += _r9;
        double _r10 = _t9 * _r_d15;
        _d_xSqrtTime += _r10;
        _d_xDen -= _r_d15;
    }
    {
        double _r_d14 = _d_xD1;
        _d_xD1 += _r_d14;
        _d_xLogTerm += _r_d14;
        _d_xD1 -= _r_d14;
    }
    {
        double _r_d13 = _d_xD1;
        double _r7 = _r_d13 * _t6;
        _d_xD1 += _r7;
        double _r8 = _t7 * _r_d13;
        _d_xTime += _r8;
        _d_xD1 -= _r_d13;
    }
    {
        double _r_d12 = _d_xD1;
        _d_xRiskFreeRate += _r_d12;
        _d_xPowerTerm += _r_d12;
        _d_xD1 -= _r_d12;
    }
    {
        double _r_d11 = _d_xPowerTerm;
        double _r6 = _r_d11 * 0.5;
        _d_xPowerTerm += _r6;
        _d_xPowerTerm -= _r_d11;
    }
    {
        double _r_d10 = _d_xPowerTerm;
        double _r4 = _r_d10 * _t4;
        _d_xVolatility += _r4;
        double _r5 = _t5 * _r_d10;
        _d_xVolatility += _r5;
        _d_xPowerTerm -= _r_d10;
    }
    {
        double _r_d9 = _d_xLogTerm;
        _d_logValues += _r_d9;
        _d_xLogTerm -= _r_d9;
    }
    {
        double _r_d8 = _d_logValues;
        double _r3 = _r_d8 * clad::custom_derivatives::log_pushforward(_t3, 1.).pushforward;
        _d_clad_log_logValues_ += _r3;
        _d_logValues -= _r_d8;
    }
    {
        double _r_d7 = _d_clad_log_logValues_;
        double _r1 = _r_d7 / _t1;
        * _d_sptprice += _r1;
        double _r2 = _r_d7 * -_t2 / (_t1 * _t1);
        * _d_strike += _r2;
        _delta_clad_log_logValues_ += clad::doApprox(_r_d7, _EERepl_clad_log_logValues_1, "clad_log_logValues_");
        _d_clad_log_logValues_ -= _r_d7;
    }
    {
        double _r_d6 = _d_xSqrtTime;
        double _r0 = _r_d6 * clad::custom_derivatives::sqrt_pushforward(_t0, 1.).pushforward;
        _d_clad_sqr_xSqrtTime_ += _r0;
        _d_xSqrtTime -= _r_d6;
    }
    {
        double _r_d5 = _d_clad_sqr_xSqrtTime_;
        _d_xTime += _r_d5;
        _delta_clad_sqr_xSqrtTime_ += clad::doApprox(_r_d5, _EERepl_clad_sqr_xSqrtTime_, "clad_sqr_xSqrtTime_");
        _d_clad_sqr_xSqrtTime_ -= _r_d5;
    }
    {
        double _r_d4 = _d_xTime;
        * _d_time += _r_d4;
        _d_xTime -= _r_d4;
    }
    {
        double _r_d3 = _d_xVolatility;
        * _d_volatility += _r_d3;
        _d_xVolatility -= _r_d3;
    }
    {
        double _r_d2 = _d_xRiskFreeRate;
        * _d_rate += _r_d2;
        _d_xRiskFreeRate -= _r_d2;
    }
    {
        double _r_d1 = _d_xStrikePrice;
        * _d_strike += _r_d1;
        _d_xStrikePrice -= _r_d1;
    }
    {
        double _r_d0 = _d_xStockPrice;
        * _d_sptprice += _r_d0;
        _d_xStockPrice -= _r_d0;
    }
    double _delta_sptprice = 0;
    double _delta_strike = 0;
    double _delta_rate = 0;
    double _delta_volatility = 0;
    double _delta_time = 0;
    _final_error += _delta_volatility + _delta_rate + _delta_sptprice + _delta_xLocal_3 + _delta_inv_sqrt_2xPI + _delta_xLocal + _delta_xK2_4 + _delta_xK2_2 + _delta_xK2 + _delta_expValues + _delta_xNPrimeofX + _delta_clad_exp_FutureValueX_ + _delta_xD2 + _delta_clad_exp_expValues_ + _delta_xD1 + _delta_clad_sqr_xSqrtTime_ + _delta_xLogTerm + _delta_clad_log_logValues_ + _delta_logValues + _delta_strike + _delta_xSqrtTime + _delta_xTime + _delta_xStockPrice + _delta_OptionPrice + _delta_xK2_5 + _delta_xStrikePrice + _delta_xRiskFreeRate + _delta_time + _delta_xVolatility + _delta_xK2_3 + _delta_xPowerTerm + _delta_xDen + _delta_d1 + _delta_d2 + _delta_FutureValueX + _delta_NofXd1 + _delta_xLocal_2 + _delta_NofXd2 + _delta_NegNofXd1 + _delta_NegNofXd2 + _delta_InputX + _delta_OutputX + _delta_xLocal_1 + _delta_xInput;
}
}
