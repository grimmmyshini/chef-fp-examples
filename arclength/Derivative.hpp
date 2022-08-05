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
        
        clad::tape<int> _t2 = {};
        
        clad::tape<double> _t3 = {};
        double _delta_x = 0;
        
        clad::tape<double> _EERepl_x1 = {};
        clad::tape<double> _t4 = {};
        double _delta_d2 = 0;
        
        clad::tape<double> _EERepl_d21 = {};
        clad::tape<double> _t5 = {};
        clad::tape<unsigned long> _t6 = {};
        clad::tape<int> _t7 = {};
        
        clad::tape<int> _t8 = {};
        
        clad::tape<double> _t9 = {};
        clad::tape<double> _EERepl_d22 = {};
        clad::tape<double> _t10 = {};
        double _delta_t3 = 0;
        
        clad::tape<double> _t11 = {};
        clad::tape<double> _EERepl_t31 = {};
        clad::tape<double> _t12 = {};
        clad::tape<double> _t13 = {};
        double _delta_s1 = 0;
        double _EERepl_s10 = s1;
        clad::tape<double> _t14 = {};
        clad::tape<double> _t15 = {};
        clad::tape<double> _t16 = {};
        clad::tape<double> _EERepl_s11 = {};
        clad::tape<double> _t17 = {};
        int i, k;
        double t2, h = 3.14159265358979319992L / 1000000, x, d2, t3;
        _EERepl_t24 = t2;
        _EERepl_t23 = t2;
        _EERepl_t22 = t2;
        _EERepl_t21 = t2;
        _EERepl_t20 = t2;
        _t0 = 0;
        {
            for (i = 1; i <= 1000000; clad::push(_t2, i), (i += 1))
            {
                _t0++;
                clad::push(_t3, x);
                x = i * h;
                clad::push(_EERepl_x1, x);
                clad::push(_t4, d2);
                d2 = 1.;
                clad::push(_EERepl_d21, d2);
                clad::push(_t5, t3);
                t3 = x;
                clad::push(_t6, 0UL);
                {
                    clad::push(_t7, k);
                    for (k = 1; k <= 5; clad::push(_t8, k), (k += 1))
                    {
                        clad::back(_t6)++;
                        clad::push(_t9, d2);
                        d2 = 2. * d2;
                        clad::push(_EERepl_d22, d2);
                        clad::push(_t10, t3);
                        t3 = t3 + sin(clad::push(_t11, d2 * x)) / d2;
                        clad::push(_EERepl_t31, t3);
                    }
                }
                clad::push(_t12, t2);
                t2 = t3;
                clad::push(_t13, s1);
                s1 = s1 + sqrt(clad::push(_t16, h * h + clad::push(_t15, (t2 - t1)) * clad::push(_t14, (t2 - t1))));
                clad::push(_EERepl_s11, s1);
                clad::push(_t17, t1);
                t1 = t2;
            }
        }
        double do_fun_return = s1;
        goto _label0;
    _label0:
        *_d_s1 += 1;
        for (; _t0; _t0--)
        {
            {
                i = clad::pop(_t2);
                int _r_d1 = _d_i;
                _d_i += _r_d1;
                _d_i -= _r_d1;
            }
            {
                {
                    t1 = clad::pop(_t17);
                    double _r_d11 = *_d_t1;
                    _d_t2 += _r_d11;
                    *_d_t1 -= _r_d11;
                    *_d_t1;
                }
                {
                    s1 = clad::pop(_t13);
                    double _r_d10 = *_d_s1;
                    *_d_s1 += _r_d10;
                    double _r14 = _r_d10 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t16), 1.).pushforward;
                    double _r15 = _r14 * h;
                    _d_h += _r15;
                    double _r16 = h * _r14;
                    _d_h += _r16;
                    double _r17 = _r14 * clad::pop(_t14);
                    _d_t2 += _r17;
                    *_d_t1 += -_r17;
                    double _r18 = clad::pop(_t15) * _r14;
                    _d_t2 += _r18;
                    *_d_t1 += -_r18;
                    double _r19 = clad::pop(_EERepl_s11);
                    _delta_s1 += clad::getErrorVal(_r_d10, _r19, "s1");
                    *_d_s1 -= _r_d10;
                    *_d_s1;
                }
                {
                    t2 = clad::pop(_t12);
                    double _r_d9 = _d_t2;
                    _d_t3 += _r_d9;
                    _d_t2 -= _r_d9;
                }
                {
                    for (; clad::back(_t6); clad::back(_t6)--)
                    {
                        {
                            k = clad::pop(_t8);
                            int _r_d6 = _d_k;
                            _d_k += _r_d6;
                            _d_k -= _r_d6;
                        }
                        {
                            {
                                t3 = clad::pop(_t10);
                                double _r_d8 = _d_t3;
                                _d_t3 += _r_d8;
                                double _r8 = _r_d8 / d2;
                                double _r9 = _r8 * clad::custom_derivatives::sin_pushforward(clad::pop(_t11), 1.).pushforward;
                                double _r10 = _r9 * x;
                                _d_d2 += _r10;
                                double _r11 = d2 * _r9;
                                _d_x += _r11;
                                double _r12 = _r_d8 * -sin(clad::push(_t11, d2 * x)) / (d2 * d2);
                                _d_d2 += _r12;
                                double _r13 = clad::pop(_EERepl_t31);
                                _delta_t3 += clad::getErrorVal(_r_d8, _r13, "t3");
                                _d_t3 -= _r_d8;
                            }
                            {
                                d2 = clad::pop(_t9);
                                double _r_d7 = _d_d2;
                                double _r5 = _r_d7 * d2;
                                double _r6 = 2. * _r_d7;
                                _d_d2 += _r6;
                                double _r7 = clad::pop(_EERepl_d22);
                                _delta_d2 += clad::getErrorVal(_r_d7, _r7, "d2");
                                _d_d2 -= _r_d7;
                            }
                        }
                    }
                    clad::pop(_t6);
                }
                {
                    t3 = clad::pop(_t5);
                    double _r_d4 = _d_t3;
                    _d_x += _r_d4;
                    _d_t3 -= _r_d4;
                }
                {
                    d2 = clad::pop(_t4);
                    double _r_d3 = _d_d2;
                    double _r4 = clad::pop(_EERepl_d21);
                    _delta_d2 += clad::getErrorVal(_r_d3, _r4, "d2");
                    _d_d2 -= _r_d3;
                }
                {
                    x = clad::pop(_t3);
                    double _r_d2 = _d_x;
                    double _r1 = _r_d2 * h;
                    _d_i += _r1;
                    double _r2 = i * _r_d2;
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

} // clad