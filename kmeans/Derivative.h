namespace clad {

inline void euclid_dist_2_grad(double* pt1, double* pt2, int numdims,
                               clad::array_ref<double> _d_pt1,
                               clad::array_ref<double> _d_pt2,
                               clad::array_ref<int> _d_numdims,
                               double& _final_error) {
  int _d_i = 0;
  double _d_ans = 0;
  double _delta_ans = 0;
  double _EERepl_ans0;
  double _d_arg = 0;
  double _delta_arg = 0;
  double _EERepl_arg0;
  unsigned long _t0;
  clad::tape<int> _t1 = {};
  clad::tape<int> _t3 = {};
  clad::tape<double> _EERepl_arg1 = {};
  clad::tape<double> _t5 = {};
  clad::tape<double> _t6 = {};
  clad::tape<double> _EERepl_ans1 = {};
  int i;
  double ans = 0.;
  _EERepl_ans0 = ans;
  double arg;
  _EERepl_arg0 = arg;
  _t0 = 0;
  for (i = 0; i < numdims; i++) {
    _t0++;
    arg = (pt1[clad::push(_t1, i)] - pt2[clad::push(_t3, i)]);
    clad::push(_EERepl_arg1, arg);
    ans += clad::push(_t6, arg) * clad::push(_t5, arg);
    clad::push(_EERepl_ans1, ans);
  }
  double euclid_dist_2_return = ans;
  goto _label0;
_label0:
  _d_ans += 1;
  for (; _t0; _t0--) {
    {
      double _r_d2 = _d_ans;
      _d_ans += _r_d2;
      double _r1 = _r_d2 * clad::pop(_t5);
      _d_arg += _r1;
      double _r2 = clad::pop(_t6) * _r_d2;
      _d_arg += _r2;
      double _r3 = clad::pop(_EERepl_ans1);
      _delta_ans += getErrorVal(_r_d2, _r3, "ans");
      _d_ans -= _r_d2;
    }
    {
      double _r_d1 = _d_arg;
      int _t2 = clad::pop(_t1);
      _d_pt1[_t2] += _r_d1;
      int _t4 = clad::pop(_t3);
      _d_pt2[_t4] += -_r_d1;
      double _r0 = clad::pop(_EERepl_arg1);
      _delta_arg += getErrorVal(_r_d1, _r0, "arg");
      _d_arg -= _r_d1;
    }
  }
  _delta_ans += getErrorVal(_d_ans, _EERepl_ans0, "ans");
  clad::array<double> _delta_pt1(_d_pt1.size());
  int i0 = 0;
  for (; i0 < _d_pt1.size(); i0++) {
    double _t7 = getErrorVal(_d_pt1[i0], pt1[i0], "attributes");
    _delta_pt1[i0] += _t7;
    _final_error += _t7;
  }
  clad::array<double> _delta_pt2(_d_pt2.size());
  i0 = 0;
  for (; i0 < _d_pt2.size(); i0++) {
    double _t8 = getErrorVal(_d_pt2[i0], pt2[i0], "clusters");
    _delta_pt2[i0] += _t8;
    _final_error += _t8;
  }
  _final_error += _delta_ans + _delta_arg;
}

} // namespace clad
