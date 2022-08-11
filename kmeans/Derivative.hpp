namespace clad
{

    inline void euclid_dist_2_grad(double *pt1, double *pt2, int numdims, clad::array_ref<double> _d_pt1, clad::array_ref<double> _d_pt2, clad::array_ref<int> _d_numdims, double &_final_error)
    {
        int _d_i = 0;
        double _d_ans = 0;
        double _delta_ans = 0;
        double _EERepl_ans0;
        double _d_x = 0, _d_y = 0;
        double _delta_x = 0;
        double _EERepl_x0;
        double _delta_x0 = 0;
        double _EERepl_x1;
        unsigned long _t0;

        clad::tape<int> _t1 = {};
        clad::tape<int> _t3 = {};
        clad::tape<double> _t5 = {};
        clad::tape<double> _t6 = {};
        clad::tape<double> _EERepl_ans1 = {};
        int i;
        double ans = 0.;
        _EERepl_ans0 = ans;
        double x, y;
        _EERepl_x1 = x;
        _EERepl_x0 = x;
        _t0 = 0;
        for (i = 0; i < numdims; i++)
        {
            _t0++;
            x = pt1[clad::push(_t1, i)];
            y = pt2[clad::push(_t3, i)];
            ans += clad::push(_t6, (x - y)) * clad::push(_t5, (x - y));
            clad::push(_EERepl_ans1, ans);
        }
        double euclid_dist_2_return = ans;
        goto _label0;
    _label0:
        _d_ans += 1;
        for (; _t0; _t0--)
        {
            {
                double _r_d3 = _d_ans;
                _d_ans += _r_d3;
                double _r0 = _r_d3 * clad::pop(_t5);
                _d_x += _r0;
                _d_y += -_r0;
                double _r1 = clad::pop(_t6) * _r_d3;
                _d_x += _r1;
                _d_y += -_r1;
                double _r2 = clad::pop(_EERepl_ans1);
                _delta_ans += clad::getErrorVal(_r_d3, _r2, "diff");
                _d_ans -= _r_d3;
            }
            {
                double _r_d2 = _d_y;
                int _t4 = clad::pop(_t3);
                _d_pt2[_t4] += _r_d2;
                _d_y -= _r_d2;
            }
            {
                double _r_d1 = _d_x;
                int _t2 = clad::pop(_t1);
                _d_pt1[_t2] += _r_d1;
                _d_x -= _r_d1;
            }
        }
        _delta_ans += clad::getErrorVal(_d_ans, _EERepl_ans0, "diff");
        clad::array<double> _delta_pt1(_d_pt1.size());
        i = 0;
        for (; i < _d_pt1.size(); i++)
        {
            double _t7 = clad::getErrorVal(_d_pt1[i], pt1[i], "attributes");
            _delta_pt1[i] += _t7;
            _final_error += _t7;
        }
        clad::array<double> _delta_pt2(_d_pt2.size());
        i = 0;
        for (; i < _d_pt2.size(); i++)
        {
            double _t8 = clad::getErrorVal(_d_pt2[i], pt2[i], "clusters");
            _delta_pt2[i] += _t8;
            _final_error += _t8;
        }
        _final_error += _delta_x + _delta_ans;
    }

} // namespace clad