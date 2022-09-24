#include <cmath>
#include <iostream>

namespace clad {
void fun_pullback(double x, double _d_y, clad::array_ref<double> _d_x,
                  double& _final_error) {
  double _delta_d1 = 0;
  double _EERepl_d10;
  double _delta_d2 = 0;
  double _EERepl_d20;
  double _delta_t3 = 0;
  double _EERepl_t30;
  unsigned long _t0;
  clad::tape<double> _t1 = {};
  clad::tape<double> _EERepl_d21 = {};
  clad::tape<double> _t2 = {};
  clad::tape<double> _t3 = {};
  clad::tape<double> _t4 = {};
  clad::tape<double> _t5 = {};
  clad::tape<double> _t6 = {};
  clad::tape<double> _EERepl_t31 = {};
  double _d_d1 = 0;
  double d1 = 1.;
  _EERepl_d10 = d1;
  double _d_d2 = 0;
  double d2 = d1;
  _EERepl_d20 = d2;
  double _d_t3 = 0;
  double t3 = x;
  _EERepl_t30 = t3;
  int _d_k = 0;
  int k;
  _t0 = 0;
  for (k = 1; k <= 5; k += 1) {
    _t0++;
    d2 = 2. * clad::push(_t1, d2);
    clad::push(_EERepl_d21, d2);
    t3 = t3 + clad::push(_t6, sin(clad::push(_t5, clad::push(_t4, d2) *
                                                      clad::push(_t3, x)))) /
                  clad::push(_t2, d2);
    clad::push(_EERepl_t31, t3);
  }
  double fun_return = t3;
  goto _label0;
_label0:
  _d_t3 += _d_y;
  for (; _t0; _t0--) {
    {
      int _r_d1 = _d_k;
      _d_k += _r_d1;
      _d_k -= _r_d1;
    }
    {
      {
        double _r_d3 = _d_t3;
        _d_t3 += _r_d3;
        double _r3 = clad::pop(_t2);
        double _r4 = _r_d3 / _r3;
        double _r5 =
            _r4 * clad::custom_derivatives::sin_pushforward(clad::pop(_t5), 1.)
                      .pushforward;
        double _r6 = _r5 * clad::pop(_t3);
        _d_d2 += _r6;
        double _r7 = clad::pop(_t4) * _r5;
        *_d_x += _r7;
        double _r8 = _r_d3 * -clad::pop(_t6) / (_r3 * _r3);
        _d_d2 += _r8;
        double _r9 = clad::pop(_EERepl_t31);
        _delta_t3 += clad::getErrorVal(_r_d3, _r9, "t3");
        _d_t3 -= _r_d3;
      }
      {
        double _r_d2 = _d_d2;
        double _r0 = _r_d2 * clad::pop(_t1);
        double _r1 = 2. * _r_d2;
        _d_d2 += _r1;
        double _r2 = clad::pop(_EERepl_d21);
        _delta_d2 += clad::getErrorVal(_r_d2, _r2, "d2");
        _d_d2 -= _r_d2;
      }
    }
  }
  *_d_x += _d_t3;
  _d_d1 += _d_d2;
  _delta_d1 += clad::getErrorVal(_d_d1, _EERepl_d10, "d1");
  double _delta_x = 0;
  _delta_x += clad::getErrorVal(*_d_x, x, "x");
  _final_error += _delta_t3 + _delta_x + _delta_d1 + _delta_d2;
}

void do_fun_grad(int iterations, clad::array_ref<int> _d_iterations,
                 double& _final_error) {
  double _delta_t2 = 0;
  double _EERepl_t20;
  double _t0;
  double _delta_h = 0;
  double _EERepl_h0;
  double _delta_t1 = 0;
  double _EERepl_t10;
  double _delta_s1 = 0;
  double _EERepl_s10;
  unsigned long _t1;
  clad::tape<double> _t2 = {};
  clad::tape<int> _t3 = {};
  clad::tape<double> _t4 = {};
  clad::tape<double> _EERepl_t21 = {};
  clad::tape<double> _t6 = {};
  clad::tape<double> _t7 = {};
  clad::tape<double> _t8 = {};
  clad::tape<double> _t9 = {};
  clad::tape<double> _t10 = {};
  clad::tape<double> _EERepl_s11 = {};
  clad::tape<double> _EERepl_t11 = {};
  double _d_t2 = 0;
  double t2;
  _EERepl_t20 = t2;
  _t0 = iterations;
  double _d_h = 0;
  double h = M_PI / _t0;
  _EERepl_h0 = h;
  double _d_t1 = 0;
  double t1 = 0;
  _EERepl_t10 = t1;
  double _d_s1 = 0;
  double s1 = 0;
  _EERepl_s10 = s1;
  int _d_i = 0;
  int i;
  _t1 = 0;
  for (i = 1; i <= iterations; i += 1) {
    _t1++;
    t2 = fun(clad::push(_t4, clad::push(_t3, i) * clad::push(_t2, h)));
    clad::push(_EERepl_t21, t2);
    s1 = s1 + sqrt(clad::push(_t10, clad::push(_t7, h) * clad::push(_t6, h) +
                                        clad::push(_t9, (t2 - t1)) *
                                            clad::push(_t8, (t2 - t1))));
    clad::push(_EERepl_s11, s1);
    t1 = +t2;
    clad::push(_EERepl_t11, t1);
  }
  double do_fun_return = s1;
  goto _label0;
_label0:
  _d_s1 += 1;
  for (; _t1; _t1--) {
    {
      int _r_d1 = _d_i;
      _d_i += _r_d1;
      _d_i -= _r_d1;
    }
    {
      {
        double _r_d4 = _d_t1;
        _d_t2 += _r_d4;
        double _r12 = clad::pop(_EERepl_t11);
        _delta_t1 += clad::getErrorVal(_r_d4, _r12, "t1");
        _d_t1 -= _r_d4;
      }
      {
        double _r_d3 = _d_s1;
        _d_s1 += _r_d3;
        double _r6 = _r_d3 * clad::custom_derivatives::sqrt_pushforward(
                                 clad::pop(_t10), 1.)
                                 .pushforward;
        double _r7 = _r6 * clad::pop(_t6);
        _d_h += _r7;
        double _r8 = clad::pop(_t7) * _r6;
        _d_h += _r8;
        double _r9 = _r6 * clad::pop(_t8);
        _d_t2 += _r9;
        _d_t1 += -_r9;
        double _r10 = clad::pop(_t9) * _r6;
        _d_t2 += _r10;
        _d_t1 += -_r10;
        double _r11 = clad::pop(_EERepl_s11);
        _delta_s1 += clad::getErrorVal(_r_d3, _r11, "s1");
        _d_s1 -= _r_d3;
      }
      {
        double _r_d2 = _d_t2;
        double _grad0 = 0.;
        double _t5 = 0;
        fun_pullback(clad::pop(_t4), _r_d2, &_grad0, _t5);
        double _r2 = _grad0;
        double _r3 = _r2 * clad::pop(_t2);
        _d_i += _r3;
        double _r4 = clad::pop(_t3) * _r2;
        _d_h += _r4;
        double _r5 = clad::pop(_EERepl_t21);
        _delta_t2 += _t5;
        _d_t2 -= _r_d2;
      }
    }
  }
  _delta_s1 += clad::getErrorVal(_d_s1, _EERepl_s10, "s1");
  _delta_t1 += clad::getErrorVal(_d_t1, _EERepl_t10, "t1");
  {
    double _r0 = _d_h / _t0;
    double _r1 = _d_h * -M_PI / (_t0 * _t0);
    *_d_iterations += _r1;
  }
  _final_error += _delta_s1 + _delta_t1 + _delta_t2 + _delta_h;
}
} // namespace clad
