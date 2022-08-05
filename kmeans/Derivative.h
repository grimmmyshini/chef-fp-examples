inline void euclid_dist_2_grad(float *pt1, float *pt2, int numdims, clad::array_ref<float> _d_pt1, clad::array_ref<float> _d_pt2, clad::array_ref<int> _d_numdims, double &_final_error) {
    int _d_i = 0;
    float _d_ans = 0;
    double _delta_ans = 0;
    float _EERepl_ans0;
    unsigned long _t0;
    
    clad::tape<float> _t1 = {};
    clad::tape<int> _t2 = {};
    clad::tape<int> _t4 = {};
    clad::tape<float> _t6 = {};
    clad::tape<int> _t7 = {};
    clad::tape<int> _t9 = {};
    clad::tape<float> _EERepl_ans1 = {};
    int i;
    float ans = 0.;
    _EERepl_ans0 = ans;
    _t0 = 0;
    for (i = 0; i < numdims; i++) {
        _t0++;
        ans += clad::push(_t6, (pt1[clad::push(_t2, i)] - pt2[clad::push(_t4, i)])) * clad::push(_t1, (pt1[clad::push(_t7, i)] - pt2[clad::push(_t9, i)]));
        clad::push(_EERepl_ans1, ans);
    }
    float euclid_dist_2_return = ans;
    goto _label0;
  _label0:
    _d_ans += 1;
    for (; _t0; _t0--) {
        float _r_d1 = _d_ans;
        _d_ans += _r_d1;
        float _r0 = _r_d1 * clad::pop(_t1);
        int _t3 = clad::pop(_t2);
        _d_pt1[_t3] += _r0;
        int _t5 = clad::pop(_t4);
        _d_pt2[_t5] += -_r0;
        float _r1 = clad::pop(_t6) * _r_d1;
        int _t8 = clad::pop(_t7);
        _d_pt1[_t8] += _r1;
        int _t10 = clad::pop(_t9);
        _d_pt2[_t10] += -_r1;
        float _r2 = clad::pop(_EERepl_ans1);
        _delta_ans += clad::getErrorVal(_r_d1, _r2, "ans");
        _d_ans -= _r_d1;
    }
    _delta_ans += clad::getErrorVal(_d_ans, _EERepl_ans0, "ans");
    clad::array<float> _delta_pt1(_d_pt1.size());
    i = 0;
    for (; i < _d_pt1.size(); i++) {
        double _t11 = clad::getErrorVal(_d_pt1[i], pt1[i], "pt1");
        _delta_pt1[i] += _t11;
        _final_error += _t11;
    }
    clad::array<float> _delta_pt2(_d_pt2.size());
    i = 0;
    for (; i < _d_pt2.size(); i++) {
        double _t12 = clad::getErrorVal(_d_pt2[i], pt2[i], "pt2");
        _delta_pt2[i] += _t12;
        _final_error += _t12;
    }
    _final_error += _delta_ans;
}
