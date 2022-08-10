namespace clad
{
    void do_fun_grad(double s1, double t1, clad::array_ref<double> _d_s1, clad::array_ref<double> _d_t1, double &_final_error)
    {
        int _d_i = 0, _d_k = 0;
        double _d_t2 = 0, _d_h = 0, _d_x = 0, _d_d2 = 0, _d_t3 = 0;
        double _delta_t2 = 0;
        double _EERepl_t20;
        double _delta_t20 = 0;
        double _EERepl_t21;
        double _delta_t21 = 0;
        double _EERepl_t22;
        double _delta_t22 = 0;
        double _EERepl_t23;
        double _delta_t23 = 0;
        double _EERepl_t24;
        unsigned long _t0;
        double _delta_x = 0;
        clad::tape<double> _t1 = {};
        clad::tape<int> _t2 = {};
        clad::tape<double> _EERepl_x1 = {};
        double _delta_d2 = 0;
        clad::tape<double> _EERepl_d21 = {};
        clad::tape<unsigned long> _t3 = {};
        clad::tape<double> _t4 = {};
        clad::tape<double> _EERepl_d22 = {};
        double _delta_t3 = 0;
        clad::tape<double> _t5 = {};
        clad::tape<double> _t6 = {};
        clad::tape<double> _t7 = {};
        clad::tape<double> _t8 = {};
        clad::tape<double> _t9 = {};
        clad::tape<double> _EERepl_t31 = {};
        double _delta_s1 = 0;
        double _EERepl_s10 = s1;
        clad::tape<double> _t10 = {};
        clad::tape<double> _t11 = {};
        clad::tape<double> _t12 = {};
        clad::tape<double> _t13 = {};
        clad::tape<double> _t14 = {};
        clad::tape<double> _EERepl_s11 = {};
        int i, k;
        double t2, h = 3.14159265358979319992L / 1000000, x, d2, t3;
        _EERepl_t24 = t2;
        _EERepl_t23 = t2;
        _EERepl_t22 = t2;
        _EERepl_t21 = t2;
        _EERepl_t20 = t2;
        _t0 = 0;
        for (i = 1; i <= 1000000; i += 1)
        {
            _t0++;
            x = clad::push(_t2, i) * clad::push(_t1, h);
            clad::push(_EERepl_x1, x);
            d2 = 1.;
            clad::push(_EERepl_d21, d2);
            t3 = x;
            clad::push(_t3, 0UL);
            for (k = 1; k <= 5; k += 1)
            {
                clad::back(_t3)++;
                d2 = 2. * clad::push(_t4, d2);
                clad::push(_EERepl_d22, d2);
                t3 = t3 + clad::push(_t9, sin(clad::push(_t8, clad::push(_t7, d2) * clad::push(_t6, x)))) / clad::push(_t5, d2);
                clad::push(_EERepl_t31, t3);
            }
            t2 = t3;
            s1 = s1 + sqrt(clad::push(_t14, clad::push(_t11, h) * clad::push(_t10, h) + clad::push(_t13, (t2 - t1)) * clad::push(_t12, (t2 - t1))));
            clad::push(_EERepl_s11, s1);
            t1 = t2;
        }
        double do_fun_return = s1;
        goto _label0;
    _label0:
        *_d_s1 += 1;
        for (; _t0; _t0--)
        {
            {
                int _r_d1 = _d_i;
                _d_i += _r_d1;
                _d_i -= _r_d1;
            }
            {
                {
                    double _r_d11 = *_d_t1;
                    _d_t2 += _r_d11;
                    *_d_t1 -= _r_d11;
                    *_d_t1;
                }
                {
                    double _r_d10 = *_d_s1;
                    *_d_s1 += _r_d10;
                    double _r15 = _r_d10 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t14), 1.).pushforward;
                    double _r16 = _r15 * clad::pop(_t10);
                    _d_h += _r16;
                    double _r17 = clad::pop(_t11) * _r15;
                    _d_h += _r17;
                    double _r18 = _r15 * clad::pop(_t12);
                    _d_t2 += _r18;
                    *_d_t1 += -_r18;
                    double _r19 = clad::pop(_t13) * _r15;
                    _d_t2 += _r19;
                    *_d_t1 += -_r19;
                    double _r20 = clad::pop(_EERepl_s11);
                    _delta_s1 += clad::getErrorVal(_r_d10, _r20, "s1");
                    *_d_s1 -= _r_d10;
                    *_d_s1;
                }
                {
                    double _r_d9 = _d_t2;
                    _d_t3 += _r_d9;
                    _d_t2 -= _r_d9;
                }
                {
                    for (; clad::back(_t3); clad::back(_t3)--)
                    {
                        {
                            int _r_d6 = _d_k;
                            _d_k += _r_d6;
                            _d_k -= _r_d6;
                        }
                        {
                            {
                                double _r_d8 = _d_t3;
                                _d_t3 += _r_d8;
                                double _r8 = clad::pop(_t5);
                                double _r9 = _r_d8 / _r8;
                                double _r10 = _r9 * clad::custom_derivatives::sin_pushforward(clad::pop(_t8), 1.).pushforward;
                                double _r11 = _r10 * clad::pop(_t6);
                                _d_d2 += _r11;
                                double _r12 = clad::pop(_t7) * _r10;
                                _d_x += _r12;
                                double _r13 = _r_d8 * -clad::pop(_t9) / (_r8 * _r8);
                                _d_d2 += _r13;
                                double _r14 = clad::pop(_EERepl_t31);
                                _delta_t3 += clad::getErrorVal(_r_d8, _r14, "t3");
                                _d_t3 -= _r_d8;
                            }
                            {
                                double _r_d7 = _d_d2;
                                double _r5 = _r_d7 * clad::pop(_t4);
                                double _r6 = 2. * _r_d7;
                                _d_d2 += _r6;
                                double _r7 = clad::pop(_EERepl_d22);
                                _delta_d2 += clad::getErrorVal(_r_d7, _r7, "d2");
                                _d_d2 -= _r_d7;
                            }
                        }
                    }
                    clad::pop(_t3);
                }
                {
                    double _r_d4 = _d_t3;
                    _d_x += _r_d4;
                    _d_t3 -= _r_d4;
                }
                {
                    double _r_d3 = _d_d2;
                    double _r4 = clad::pop(_EERepl_d21);
                    _delta_d2 += clad::getErrorVal(_r_d3, _r4, "d2");
                    _d_d2 -= _r_d3;
                }
                {
                    double _r_d2 = _d_x;
                    double _r1 = _r_d2 * clad::pop(_t1);
                    _d_i += _r1;
                    double _r2 = clad::pop(_t2) * _r_d2;
                    _d_h += _r2;
                    double _r3 = clad::pop(_EERepl_x1);
                    _delta_x += clad::getErrorVal(_r_d2, _r3, "x");
                    _d_x -= _r_d2;
                }
            }
        }
        long double _r0 = _d_h / 1000000;
        _delta_s1 += clad::getErrorVal(*_d_s1, _EERepl_s10, "s1");
        double _delta_t1 = 0;
        _delta_t1 += clad::getErrorVal(*_d_t1, t1, "t1");
        _final_error += _delta_t1 + _delta_s1 + _delta_t3 + _delta_d2 + _delta_x + _delta_t2;
    }
}