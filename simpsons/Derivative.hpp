namespace clad
{
    void simpsons_grad(double a, double b, clad::array_ref<double> _d_a, clad::array_ref<double> _d_b, double &_final_error)
    {
        int _d_n = 0;
        double _d_pi = 0;
        double _delta_pi = 0;
        double _EERepl_pi0;
        double _t0;
        double _t1;
        double _t2;
        double _d_h = 0;
        double _delta_h = 0;
        double _EERepl_h0;
        double _d_x = 0;
        double _delta_x = 0;
        double _EERepl_x0;
        double _d_tmp = 0, _d_x_pi = 0, _d_sin_x_pi = 0;
        double _delta_tmp = 0;
        double _EERepl_tmp0;
        double _delta_tmp0 = 0;
        double _EERepl_tmp1;
        double _delta_tmp1 = 0;
        double _EERepl_tmp2;
        double _delta_x_pi = 0;
        double _t3;
        double _t4;
        double _EERepl_x_pi1;
        double _delta_sin_x_pi = 0;
        double _t5;
        double _EERepl_sin_x_pi1;
        double _d_fa = 0;
        double _delta_fa = 0;
        double _EERepl_fa0;
        double _t6;
        double _t7;
        double _EERepl_x_pi2;
        double _t8;
        double _EERepl_sin_x_pi2;
        double _d_fb = 0;
        double _delta_fb = 0;
        double _EERepl_fb0;
        double _d_s1 = 0;
        double _delta_s1 = 0;
        double _EERepl_s10;
        unsigned long _t9;
        int _d_l = 0;
        clad::tape<double> _EERepl_x1 = {};
        clad::tape<double> _t10 = {};
        clad::tape<double> _t11 = {};
        clad::tape<double> _EERepl_x_pi3 = {};
        clad::tape<double> _t12 = {};
        clad::tape<double> _EERepl_sin_x_pi3 = {};
        clad::tape<double> _t13 = {};
        clad::tape<double> _EERepl_s11 = {};
        clad::tape<double> _EERepl_x2 = {};
        clad::tape<double> _t14 = {};
        clad::tape<double> _t15 = {};
        clad::tape<double> _EERepl_x_pi4 = {};
        clad::tape<double> _t16 = {};
        clad::tape<double> _EERepl_sin_x_pi4 = {};
        clad::tape<double> _t17 = {};
        clad::tape<double> _EERepl_s12 = {};
        double _EERepl_s13;
        double _EERepl_tmp3;
        double _t18;
        double _t19;
        double _EERepl_s14;
        int n = ITERATIONS;
        double pi = 3.1415926535897931;
        _EERepl_pi0 = pi;
        _t1 = (b - a);
        _t2 = n;
        _t0 = (2. * _t2);
        double h = _t1 / _t0;
        _EERepl_h0 = h;
        double x = a;
        _EERepl_x0 = x;
        double tmp, x_pi, sin_x_pi;
        _EERepl_tmp2 = tmp;
        _EERepl_tmp1 = tmp;
        _EERepl_tmp0 = tmp;
        _t4 = a;
        _t3 = pi;
        x_pi = _t4 * _t3;
        _EERepl_x_pi1 = x_pi;
        _t5 = x_pi;
        sin_x_pi = sin(_t5);
        _EERepl_sin_x_pi1 = sin_x_pi;
        double fa = sin_x_pi;
        _EERepl_fa0 = fa;
        _t7 = b;
        _t6 = pi;
        x_pi = _t7 * _t6;
        _EERepl_x_pi2 = x_pi;
        _t8 = x_pi;
        sin_x_pi = sin(_t8);
        _EERepl_sin_x_pi2 = sin_x_pi;
        double fb = sin_x_pi;
        _EERepl_fb0 = fb;
        double s1 = fa + fb;
        _EERepl_s10 = s1;
        _t9 = 0;
        for (int l = 0; l < n; l++)
        {
            _t9++;
            x = x + h;
            clad::push(_EERepl_x1, x);
            x_pi = clad::push(_t11, x) * clad::push(_t10, pi);
            clad::push(_EERepl_x_pi3, x_pi);
            sin_x_pi = sin(clad::push(_t12, x_pi));
            clad::push(_EERepl_sin_x_pi3, sin_x_pi);
            s1 = s1 + 4. * clad::push(_t13, sin_x_pi);
            clad::push(_EERepl_s11, s1);
            x = x + h;
            clad::push(_EERepl_x2, x);
            x_pi = clad::push(_t15, x) * clad::push(_t14, pi);
            clad::push(_EERepl_x_pi4, x_pi);
            sin_x_pi = sin(clad::push(_t16, x_pi));
            clad::push(_EERepl_sin_x_pi4, sin_x_pi);
            s1 = s1 + 2. * clad::push(_t17, sin_x_pi);
            clad::push(_EERepl_s12, s1);
        }
        s1 = s1 - fb;
        _EERepl_s13 = s1;
        tmp = h * 3.1415926535897931 / 3.;
        _EERepl_tmp3 = tmp;
        _t19 = s1;
        _t18 = tmp;
        s1 = _t19 * _t18;
        _EERepl_s14 = s1;
        double simpsons_return = s1;
        goto _label0;
    _label0:
        _d_s1 += 1;
        {
            double _r_d14 = _d_s1;
            double _r30 = _r_d14 * _t18;
            _d_s1 += _r30;
            double _r31 = _t19 * _r_d14;
            _d_tmp += _r31;
            _delta_s1 += clad::getErrorVal(_r_d14, _EERepl_s14, "s1");
            _d_s1 -= _r_d14;
        }
        {
            double _r_d13 = _d_tmp;
            double _r28 = _r_d13 / 3.;
            double _r29 = _r28 * 3.1415926535897931;
            _d_h += _r29;
            _delta_tmp += clad::getErrorVal(_r_d13, _EERepl_tmp3, "tmp");
            _d_tmp -= _r_d13;
        }
        {
            double _r_d12 = _d_s1;
            _d_s1 += _r_d12;
            _d_fb += -_r_d12;
            _delta_s1 += clad::getErrorVal(_r_d12, _EERepl_s13, "s1");
            _d_s1 -= _r_d12;
        }
        for (; _t9; _t9--)
        {
            {
                double _r_d11 = _d_s1;
                _d_s1 += _r_d11;
                double _r25 = _r_d11 * clad::pop(_t17);
                double _r26 = 2. * _r_d11;
                _d_sin_x_pi += _r26;
                double _r27 = clad::pop(_EERepl_s12);
                _delta_s1 += clad::getErrorVal(_r_d11, _r27, "s1");
                _d_s1 -= _r_d11;
            }
            {
                double _r_d10 = _d_sin_x_pi;
                double _r23 = _r_d10 * clad::custom_derivatives::sin_pushforward(clad::pop(_t16), 1.).pushforward;
                _d_x_pi += _r23;
                double _r24 = clad::pop(_EERepl_sin_x_pi4);
                _delta_sin_x_pi += clad::getErrorVal(_r_d10, _r24, "sin_x_pi");
                _d_sin_x_pi -= _r_d10;
            }
            {
                double _r_d9 = _d_x_pi;
                double _r20 = _r_d9 * clad::pop(_t14);
                _d_x += _r20;
                double _r21 = clad::pop(_t15) * _r_d9;
                _d_pi += _r21;
                double _r22 = clad::pop(_EERepl_x_pi4);
                _delta_x_pi += clad::getErrorVal(_r_d9, _r22, "x_pi");
                _d_x_pi -= _r_d9;
            }
            {
                double _r_d8 = _d_x;
                _d_x += _r_d8;
                _d_h += _r_d8;
                double _r19 = clad::pop(_EERepl_x2);
                _delta_x += clad::getErrorVal(_r_d8, _r19, "x");
                _d_x -= _r_d8;
            }
            {
                double _r_d7 = _d_s1;
                _d_s1 += _r_d7;
                double _r16 = _r_d7 * clad::pop(_t13);
                double _r17 = 4. * _r_d7;
                _d_sin_x_pi += _r17;
                double _r18 = clad::pop(_EERepl_s11);
                _delta_s1 += clad::getErrorVal(_r_d7, _r18, "s1");
                _d_s1 -= _r_d7;
            }
            {
                double _r_d6 = _d_sin_x_pi;
                double _r14 = _r_d6 * clad::custom_derivatives::sin_pushforward(clad::pop(_t12), 1.).pushforward;
                _d_x_pi += _r14;
                double _r15 = clad::pop(_EERepl_sin_x_pi3);
                _delta_sin_x_pi += clad::getErrorVal(_r_d6, _r15, "sin_x_pi");
                _d_sin_x_pi -= _r_d6;
            }
            {
                double _r_d5 = _d_x_pi;
                double _r11 = _r_d5 * clad::pop(_t10);
                _d_x += _r11;
                double _r12 = clad::pop(_t11) * _r_d5;
                _d_pi += _r12;
                double _r13 = clad::pop(_EERepl_x_pi3);
                _delta_x_pi += clad::getErrorVal(_r_d5, _r13, "x_pi");
                _d_x_pi -= _r_d5;
            }
            {
                double _r_d4 = _d_x;
                _d_x += _r_d4;
                _d_h += _r_d4;
                double _r10 = clad::pop(_EERepl_x1);
                _delta_x += clad::getErrorVal(_r_d4, _r10, "x");
                _d_x -= _r_d4;
            }
        }
        {
            _d_fa += _d_s1;
            _d_fb += _d_s1;
            _delta_s1 += clad::getErrorVal(_d_s1, _EERepl_s10, "s1");
        }
        _d_sin_x_pi += _d_fb;
        {
            double _r_d3 = _d_sin_x_pi;
            double _r9 = _r_d3 * clad::custom_derivatives::sin_pushforward(_t8, 1.).pushforward;
            _d_x_pi += _r9;
            _delta_sin_x_pi += clad::getErrorVal(_r_d3, _EERepl_sin_x_pi2, "sin_x_pi");
            _d_sin_x_pi -= _r_d3;
        }
        {
            double _r_d2 = _d_x_pi;
            double _r7 = _r_d2 * _t6;
            *_d_b += _r7;
            double _r8 = _t7 * _r_d2;
            _d_pi += _r8;
            _delta_x_pi += clad::getErrorVal(_r_d2, _EERepl_x_pi2, "x_pi");
            _d_x_pi -= _r_d2;
        }
        _d_sin_x_pi += _d_fa;
        {
            double _r_d1 = _d_sin_x_pi;
            double _r6 = _r_d1 * clad::custom_derivatives::sin_pushforward(_t5, 1.).pushforward;
            _d_x_pi += _r6;
            _delta_sin_x_pi += clad::getErrorVal(_r_d1, _EERepl_sin_x_pi1, "sin_x_pi");
            _d_sin_x_pi -= _r_d1;
        }
        {
            double _r_d0 = _d_x_pi;
            double _r4 = _r_d0 * _t3;
            *_d_a += _r4;
            double _r5 = _t4 * _r_d0;
            _d_pi += _r5;
            _delta_x_pi += clad::getErrorVal(_r_d0, _EERepl_x_pi1, "x_pi");
            _d_x_pi -= _r_d0;
        }
        *_d_a += _d_x;
        {
            double _r0 = _d_h / _t0;
            *_d_b += _r0;
            *_d_a += -_r0;
            double _r1 = _d_h * -_t1 / (_t0 * _t0);
            double _r2 = _r1 * _t2;
            double _r3 = 2. * _r1;
            _d_n += _r3;
            _delta_h += clad::getErrorVal(_d_h, _EERepl_h0, "h");
        }
        _delta_pi += clad::getErrorVal(_d_pi, _EERepl_pi0, "pi");
        double _delta_a = 0;
        _delta_a += clad::getErrorVal(*_d_a, a, "a");
        double _delta_b = 0;
        _delta_b += clad::getErrorVal(*_d_b, b, "b");
        _final_error += _delta_a + _delta_s1 + _delta_fa + _delta_sin_x_pi + _delta_x_pi + _delta_b + _delta_tmp + _delta_x + _delta_fb + _delta_h + _delta_pi;
    }
}
