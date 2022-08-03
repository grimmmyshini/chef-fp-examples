void ddot_pullback(int n, double *x, double *y, double _d_y0, clad::array_ref<int> _d_n, clad::array_ref<double> _d_x, clad::array_ref<double> _d_y, double &_final_error) {
    double _d_local_result = 0;
    double _delta_local_result = 0;
    double _EERepl_local_result0;
    bool _cond0;
    unsigned long _t0;
    int _d_i = 0;
    clad::tape<double> _t1 = {};
    clad::tape<int> _t2 = {};
    clad::tape<double> _t4 = {};
    clad::tape<int> _t5 = {};
    clad::tape<double> _EERepl_local_result1 = {};
    unsigned long _t7;
    int _d_i = 0;
    clad::tape<double> _t8 = {};
    clad::tape<int> _t9 = {};
    clad::tape<double> _t11 = {};
    clad::tape<int> _t12 = {};
    clad::tape<double> _EERepl_local_result2 = {};
    double local_result = 0.;
    _EERepl_local_result0 = local_result;
    _cond0 = y == x;
    if (_cond0) {
        _t0 = 0;
        for (int i = 0; i < n; i++) {
            _t0++;
            local_result += clad::push(_t4, x[clad::push(_t2, i)]) * clad::push(_t1, x[clad::push(_t5, i)]);
            clad::push(_EERepl_local_result1, local_result);
        }
    } else {
        _t7 = 0;
        for (int i = 0; i < n; i++) {
            _t7++;
            local_result += clad::push(_t11, x[clad::push(_t9, i)]) * clad::push(_t8, y[clad::push(_t12, i)]);
            clad::push(_EERepl_local_result2, local_result);
        }
    }
    double ddot_return = local_result;
    goto _label0;
  _label0:
    _d_local_result += _d_y0;
    if (_cond0)
        for (; _t0; _t0--) {
            double _r_d0 = _d_local_result;
            _d_local_result += _r_d0;
            double _r0 = _r_d0 * clad::pop(_t1);
            int _t3 = clad::pop(_t2);
            _d_x[_t3] += _r0;
            double _r1 = clad::pop(_t4) * _r_d0;
            int _t6 = clad::pop(_t5);
            _d_x[_t6] += _r1;
            double _r2 = clad::pop(_EERepl_local_result1);
            _delta_local_result += clad::getErrorVal(_r_d0, _r2, "local_result");
            _d_local_result -= _r_d0;
        }
    else
        for (; _t7; _t7--) {
            double _r_d1 = _d_local_result;
            _d_local_result += _r_d1;
            double _r3 = _r_d1 * clad::pop(_t8);
            int _t10 = clad::pop(_t9);
            _d_x[_t10] += _r3;
            double _r4 = clad::pop(_t11) * _r_d1;
            int _t13 = clad::pop(_t12);
            _d_y[_t13] += _r4;
            double _r5 = clad::pop(_EERepl_local_result2);
            _delta_local_result += clad::getErrorVal(_r_d1, _r5, "local_result");
            _d_local_result -= _r_d1;
        }
    _delta_local_result += clad::getErrorVal(_d_local_result, _EERepl_local_result0, "local_result");
    clad::array<double> _delta_x(_d_x.size());
    int i = 0;
    for (; i < _d_x.size(); i++) {
        double _t14 = clad::getErrorVal(_d_x[i], x[i], "x");
        _delta_x[i] += _t14;
        _final_error += _t14;
    }
    clad::array<double> _delta_y(_d_y.size());
    i = 0;
    for (; i < _d_y.size(); i++) {
        double _t15 = clad::getErrorVal(_d_y[i], y[i], "y");
        _delta_y[i] += _t15;
        _final_error += _t15;
    }
    _final_error += _delta_local_result;
}

void HPCCG_residual_grad(double *b, double *x, double *xexact, double *r, double *p, double *Ap, clad::array_ref<double> _d_b, clad::array_ref<double> _d_x, clad::array_ref<double> _d_xexact, clad::array_ref<double> _d_r, clad::array_ref<double> _d_p, clad::array_ref<double> _d_Ap, double &_final_error) {
    int _d_niters = 0;
    double _d_normr = 0;
    double _delta_normr = 0;
    double _EERepl_normr0;
    int _d_max_iter = 0;
    double _d_tolerance = 0;
    double _delta_tolerance = 0;
    double _EERepl_tolerance0;
    int _d_cur_nnz = 0;
    int _d_nrow = 0;
    int _d_ncol = 0;
    double _EERepl_normr1;
    double _d_rtrans = 0;
    double _delta_rtrans = 0;
    double _EERepl_rtrans0;
    double _d_oldrtrans = 0;
    double _delta_oldrtrans = 0;
    double _EERepl_oldrtrans0;
    double _d_beta = 0;
    double _delta_beta = 0;
    double _EERepl_beta0;
    unsigned long _t0;
    int _d_i = 0;
    clad::tape<int> _t1 = {};
    clad::array<double> _delta_p(_d_p.size());
    clad::array<double> _EERepl_p0(_d_p.size());
    for (int i0 = 0; i0 < _d_p.size(); i0++) {
        _EERepl_p0[i0] = p[i0];
    }
    clad::tape<int> _t3 = {};
    clad::tape<double> _t5 = {};
    clad::tape<double> _t6 = {};
    clad::tape<int> _t7 = {};
    clad::tape<double> _EERepl_p1 = {};
    unsigned long _t9;
    int _d_i = 0;
    int _EERepl_cur_nnz0 = cur_nnz;
    clad::tape<int> _t10 = {};
    clad::tape<int> _t12 = {};
    clad::array<double> _delta_Ap(_d_Ap.size());
    clad::array<double> _EERepl_Ap0(_d_Ap.size());
    for (int i0 = 0; i0 < _d_Ap.size(); i0++) {
        _EERepl_Ap0[i0] = Ap[i0];
    }
    clad::tape<double> _EERepl_Ap1 = {};
    clad::tape<unsigned long> _t14 = {};
    int _d_j = 0;
    clad::tape<int> _t15 = {};
    clad::tape<double> _t17 = {};
    clad::tape<int> _t18 = {};
    clad::tape<int> _t20 = {};
    clad::tape<double> _t22 = {};
    clad::tape<int> _t23 = {};
    clad::tape<int> _t25 = {};
    clad::tape<int> _t27 = {};
    clad::tape<double> _EERepl_Ap2 = {};
    double _EERepl_beta1;
    unsigned long _t29;
    int _d_i = 0;
    clad::tape<int> _t30 = {};
    clad::array<double> _delta_r(_d_r.size());
    clad::array<double> _EERepl_r0(_d_r.size());
    for (int i0 = 0; i0 < _d_r.size(); i0++) {
        _EERepl_r0[i0] = r[i0];
    }
    clad::tape<int> _t32 = {};
    clad::tape<double> _t34 = {};
    clad::tape<double> _t35 = {};
    clad::tape<int> _t36 = {};
    clad::tape<double> _EERepl_r1 = {};
    int _t38;
    double *_t39;
    double *_t40;
    double _EERepl_rtrans1;
    double _t42;
    double _EERepl_normr2;
    unsigned long _t43;
    int _d_k = 0;
    clad::tape<bool> _t45 = {};
    clad::tape<double> _EERepl_beta2 = {};
    clad::tape<unsigned long> _t46 = {};
    int _d_i = 0;
    clad::tape<int> _t47 = {};
    clad::tape<int> _t49 = {};
    clad::tape<double> _t51 = {};
    clad::tape<double> _t52 = {};
    clad::tape<int> _t53 = {};
    clad::tape<double> _EERepl_p2 = {};
    clad::tape<int> _t55 = {};
    clad::tape<double *> _t56 = {};
    clad::tape<double *> _t57 = {};
    clad::tape<double> _EERepl_rtrans2 = {};
    clad::tape<double> _t59 = {};
    clad::tape<double> _t60 = {};
    clad::tape<double> _EERepl_beta3 = {};
    clad::tape<unsigned long> _t61 = {};
    int _d_i = 0;
    clad::tape<int> _t62 = {};
    clad::tape<int> _t64 = {};
    clad::tape<double> _t66 = {};
    clad::tape<double> _t67 = {};
    clad::tape<int> _t68 = {};
    clad::tape<double> _EERepl_p3 = {};
    clad::tape<double> _t70 = {};
    clad::tape<double> _EERepl_normr3 = {};
    clad::tape<unsigned long> _t71 = {};
    int _d_i = 0;
    int _EERepl_cur_nnz1 = cur_nnz;
    clad::tape<int> _t72 = {};
    clad::tape<int> _t74 = {};
    clad::tape<double> _EERepl_Ap3 = {};
    clad::tape<unsigned long> _t76 = {};
    int _d_j = 0;
    clad::tape<int> _t77 = {};
    clad::tape<double> _t79 = {};
    clad::tape<int> _t80 = {};
    clad::tape<int> _t82 = {};
    clad::tape<double> _t84 = {};
    clad::tape<int> _t85 = {};
    clad::tape<int> _t87 = {};
    clad::tape<int> _t89 = {};
    clad::tape<double> _EERepl_Ap4 = {};
    clad::tape<int> _t91 = {};
    clad::tape<double *> _t92 = {};
    clad::tape<double *> _t93 = {};
    clad::tape<double> _EERepl_beta4 = {};
    clad::tape<double> _t95 = {};
    clad::tape<double> _t96 = {};
    clad::tape<double> _EERepl_beta5 = {};
    clad::tape<unsigned long> _t97 = {};
    int _d_i = 0;
    clad::tape<int> _t98 = {};
    clad::array<double> _delta_x(_d_x.size());
    clad::array<double> _EERepl_x0(_d_x.size());
    for (int i0 = 0; i0 < _d_x.size(); i0++) {
        _EERepl_x0[i0] = x[i0];
    }
    clad::tape<int> _t100 = {};
    clad::tape<double> _t102 = {};
    clad::tape<double> _t103 = {};
    clad::tape<int> _t104 = {};
    clad::tape<double> _EERepl_x1 = {};
    clad::tape<double> _EERepl_beta6 = {};
    clad::tape<unsigned long> _t106 = {};
    int _d_i = 0;
    clad::tape<int> _t107 = {};
    clad::tape<int> _t109 = {};
    clad::tape<double> _t111 = {};
    clad::tape<double> _t112 = {};
    clad::tape<int> _t113 = {};
    clad::tape<double> _EERepl_r2 = {};
    double _d_residual = 0;
    double _delta_residual = 0;
    double _EERepl_residual0;
    double _d_diff = 0;
    double _delta_diff = 0;
    double _EERepl_diff0;
    unsigned long _t115;
    int _d_i = 0;
    clad::tape<int> _t116 = {};
    clad::tape<int> _t118 = {};
    clad::tape<double> _t120 = {};
    clad::tape<double> _EERepl_diff1 = {};
    clad::tape<bool> _t123 = {};
    int niters = 0;
    double normr = 0.;
    _EERepl_normr0 = normr;
    int max_iter = 100;
    double tolerance = 0.;
    _EERepl_tolerance0 = tolerance;
    int cur_nnz;
    int nrow = A.local_nrow;
    int ncol = A.local_ncol;
    normr = 0.;
    _EERepl_normr1 = normr;
    double rtrans = 0.;
    _EERepl_rtrans0 = rtrans;
    double oldrtrans = 0.;
    _EERepl_oldrtrans0 = oldrtrans;
    double beta = 0.;
    _EERepl_beta0 = beta;
    _t0 = 0;
    for (int i = 0; i < nrow; i++) {
        _t0++;
        p[clad::push(_t1, i)] = x[clad::push(_t3, i)] + clad::push(_t6, beta) * clad::push(_t5, x[clad::push(_t7, i)]);
        clad::push(_EERepl_p1, p[clad::push(_t1, i)]);
    }
    _t9 = 0;
    for (int i = 0; i < nrow; i++) {
        _t9++;
        cur_nnz = A.nnz_in_row[clad::push(_t10, i)];
        Ap[clad::push(_t12, i)] = 0;
        clad::push(_EERepl_Ap1, Ap[clad::push(_t12, i)]);
        clad::push(_t14, 0UL);
        for (int j = 0; j < cur_nnz; j++) {
            clad::back(_t14)++;
            Ap[clad::push(_t15, i)] += clad::push(_t22, A.ptr_to_vals_in_row[clad::push(_t18, i)][clad::push(_t20, j)]) * clad::push(_t17, p[clad::push(_t27, A.ptr_to_inds_in_row[clad::push(_t23, i)][clad::push(_t25, j)])]);
            clad::push(_EERepl_Ap2, Ap[clad::push(_t15, i)]);
        }
    }
    beta = -1.;
    _EERepl_beta1 = beta;
    _t29 = 0;
    for (int i = 0; i < nrow; i++) {
        _t29++;
        r[clad::push(_t30, i)] = b[clad::push(_t32, i)] + clad::push(_t35, beta) * clad::push(_t34, Ap[clad::push(_t36, i)]);
        clad::push(_EERepl_r1, r[clad::push(_t30, i)]);
    }
    _t38 = nrow;
    _t39 = r;
    _t40 = r;
    rtrans = ddot(_t38, r, r);
    _EERepl_rtrans1 = rtrans;
    _t42 = rtrans;
    normr = sqrt(_t42);
    _EERepl_normr2 = normr;
    _t43 = 0;
    for (int k = 1; k < max_iter && normr > tolerance; k++) {
        _t43++;
        bool _t44 = k == 1;
        {
            if (_t44) {
                beta = 0.;
                clad::push(_EERepl_beta2, beta);
                clad::push(_t46, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t46)++;
                    p[clad::push(_t47, i)] = r[clad::push(_t49, i)] + clad::push(_t52, beta) * clad::push(_t51, r[clad::push(_t53, i)]);
                    clad::push(_EERepl_p2, p[clad::push(_t47, i)]);
                }
            } else {
                oldrtrans = rtrans;
                clad::push(_t56, r);
                clad::push(_t57, r);
                rtrans = ddot(clad::push(_t55, nrow), r, r);
                clad::push(_EERepl_rtrans2, rtrans);
                beta = clad::push(_t60, rtrans) / clad::push(_t59, oldrtrans);
                clad::push(_EERepl_beta3, beta);
                clad::push(_t61, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t61)++;
                    p[clad::push(_t62, i)] = r[clad::push(_t64, i)] + clad::push(_t67, beta) * clad::push(_t66, p[clad::push(_t68, i)]);
                    clad::push(_EERepl_p3, p[clad::push(_t62, i)]);
                }
            }
            clad::push(_t45, _t44);
        }
        normr = sqrt(clad::push(_t70, rtrans));
        clad::push(_EERepl_normr3, normr);
        clad::push(_t71, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t71)++;
            cur_nnz = A.nnz_in_row[clad::push(_t72, i)];
            Ap[clad::push(_t74, i)] = 0;
            clad::push(_EERepl_Ap3, Ap[clad::push(_t74, i)]);
            clad::push(_t76, 0UL);
            for (int j = 0; j < cur_nnz; j++) {
                clad::back(_t76)++;
                Ap[clad::push(_t77, i)] += clad::push(_t84, A.ptr_to_vals_in_row[clad::push(_t80, i)][clad::push(_t82, j)]) * clad::push(_t79, p[clad::push(_t89, A.ptr_to_inds_in_row[clad::push(_t85, i)][clad::push(_t87, j)])]);
                clad::push(_EERepl_Ap4, Ap[clad::push(_t77, i)]);
            }
        }
        clad::push(_t92, p);
        clad::push(_t93, Ap);
        beta = ddot(clad::push(_t91, nrow), p, Ap);
        clad::push(_EERepl_beta4, beta);
        beta = clad::push(_t96, rtrans) / clad::push(_t95, beta);
        clad::push(_EERepl_beta5, beta);
        clad::push(_t97, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t97)++;
            x[clad::push(_t98, i)] = x[clad::push(_t100, i)] + clad::push(_t103, beta) * clad::push(_t102, p[clad::push(_t104, i)]);
            clad::push(_EERepl_x1, x[clad::push(_t98, i)]);
        }
        beta = -beta;
        clad::push(_EERepl_beta6, beta);
        clad::push(_t106, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t106)++;
            r[clad::push(_t107, i)] = r[clad::push(_t109, i)] + clad::push(_t112, beta) * clad::push(_t111, Ap[clad::push(_t113, i)]);
            clad::push(_EERepl_r2, r[clad::push(_t107, i)]);
        }
        niters = k;
    }
    double residual = 0.;
    _EERepl_residual0 = residual;
    double diff = 0;
    _EERepl_diff0 = diff;
    _t115 = 0;
    for (int i = 0; i < nrow; i++) {
        _t115++;
        diff = fabs(clad::push(_t120, x[clad::push(_t116, i)] - xexact[clad::push(_t118, i)]));
        clad::push(_EERepl_diff1, diff);
        bool _t122 = diff > residual;
        {
            if (_t122)
                residual = diff;
            clad::push(_t123, _t122);
        }
    }
    double HPCCG_residual_return = residual;
    goto _label0;
  _label0:
    _d_residual += 1;
    for (; _t115; _t115--) {
        if (clad::pop(_t123)) {
            double _r_d26 = _d_residual;
            _d_diff += _r_d26;
            _d_residual -= _r_d26;
        }
        {
            double _r_d25 = _d_diff;
            double _r64 = _r_d25 * numerical_diff::forward_central_difference(fabs, clad::pop(_t120), 0, 0, clad::pop(_t120));
            int _t117 = clad::pop(_t116);
            _d_x[_t117] += _r64;
            int _t119 = clad::pop(_t118);
            _d_xexact[_t119] += -_r64;
            double _r65 = clad::pop(_EERepl_diff1);
            _delta_diff += _t121;
            _d_diff -= _r_d25;
        }
    }
    _delta_diff += clad::getErrorVal(_d_diff, _EERepl_diff0, "diff");
    _delta_residual += clad::getErrorVal(_d_residual, _EERepl_residual0, "residual");
    for (; _t43; _t43--) {
        {
            int _r_d24 = _d_niters;
            _d_k += _r_d24;
            _d_niters -= _r_d24;
        }
        {
            for (; clad::back(_t106); clad::back(_t106)--) {
                int _t108 = clad::pop(_t107);
                double _r_d23 = _d_r[_t108];
                int _t110 = clad::pop(_t109);
                _d_r[_t110] += _r_d23;
                double _r60 = _r_d23 * clad::pop(_t111);
                _d_beta += _r60;
                double _r61 = clad::pop(_t112) * _r_d23;
                int _t114 = clad::pop(_t113);
                _d_Ap[_t114] += _r61;
                double _r62 = clad::pop(_EERepl_r2);
                int _r63 = clad::pop(_t107);
                _delta_r[_r63] += clad::getErrorVal(_r_d23, _r62, "r");
                _final_error += _delta_r[_r63];
                _d_r[_t108] -= _r_d23;
                _d_r[_t108];
            }
            clad::pop(_t106);
        }
        {
            double _r_d22 = _d_beta;
            _d_beta += -_r_d22;
            double _r59 = clad::pop(_EERepl_beta6);
            _delta_beta += clad::getErrorVal(_r_d22, _r59, "beta");
            _d_beta -= _r_d22;
        }
        {
            for (; clad::back(_t97); clad::back(_t97)--) {
                int _t99 = clad::pop(_t98);
                double _r_d21 = _d_x[_t99];
                int _t101 = clad::pop(_t100);
                _d_x[_t101] += _r_d21;
                double _r55 = _r_d21 * clad::pop(_t102);
                _d_beta += _r55;
                double _r56 = clad::pop(_t103) * _r_d21;
                int _t105 = clad::pop(_t104);
                _d_p[_t105] += _r56;
                double _r57 = clad::pop(_EERepl_x1);
                int _r58 = clad::pop(_t98);
                _delta_x[_r58] += clad::getErrorVal(_r_d21, _r57, "x");
                _final_error += _delta_x[_r58];
                _d_x[_t99] -= _r_d21;
                _d_x[_t99];
            }
            clad::pop(_t97);
        }
        {
            double _r_d20 = _d_beta;
            double _r51 = clad::pop(_t95);
            double _r52 = _r_d20 / _r51;
            _d_rtrans += _r52;
            double _r53 = _r_d20 * -clad::pop(_t96) / (_r51 * _r51);
            _d_beta += _r53;
            double _r54 = clad::pop(_EERepl_beta5);
            _delta_beta += clad::getErrorVal(_r_d20, _r54, "beta");
            _d_beta -= _r_d20;
        }
        {
            int _t88 = clad::pop(_t87);
            int _t86 = clad::pop(_t85);
            double _r_d19 = _d_beta;
            double *_r47 = clad::pop(_t92);
            double *_r49 = clad::pop(_t93);
            int _grad6 = 0;
            double _t94 = 0;
            ddot_pullback(clad::pop(_t91), _r47, _r49, _r_d19, &_grad6, _d_p, _d_Ap, _t94);
            int _r45 = _grad6;
            _d_nrow += _r45;
            clad::array<double> _r46(_d_p);
            clad::array<double> _r48(_d_Ap);
            double _r50 = clad::pop(_EERepl_beta4);
            _delta_beta += _t94;
            _d_beta -= _r_d19;
        }
        {
            for (; clad::back(_t71); clad::back(_t71)--) {
                {
                    for (; clad::back(_t76); clad::back(_t76)--) {
                        int _t78 = clad::pop(_t77);
                        double _r_d18 = _d_Ap[_t78];
                        _d_Ap[_t78] += _r_d18;
                        double _r41 = _r_d18 * clad::pop(_t79);
                        int _t81 = clad::pop(_t80);
                        int _t83 = clad::pop(_t82);
                        double _r42 = clad::pop(_t84) * _r_d18;
                        int _t90 = clad::pop(_t89);
                        _d_p[_t90] += _r42;
                        double _r43 = clad::pop(_EERepl_Ap4);
                        int _r44 = clad::pop(_t77);
                        _delta_Ap[_r44] += clad::getErrorVal(_r_d18, _r43, "Ap");
                        _final_error += _delta_Ap[_r44];
                        _d_Ap[_t78] -= _r_d18;
                        _d_Ap[_t78];
                    }
                    clad::pop(_t76);
                }
                {
                    int _t75 = clad::pop(_t74);
                    double _r_d17 = _d_Ap[_t75];
                    double _r39 = clad::pop(_EERepl_Ap3);
                    int _r40 = clad::pop(_t74);
                    _delta_Ap[_r40] += clad::getErrorVal(_r_d17, _r39, "Ap");
                    _final_error += _delta_Ap[_r40];
                    _d_Ap[_t75] -= _r_d17;
                    _d_Ap[_t75];
                }
                {
                    int _r_d16 = _d_cur_nnz;
                    int _t73 = clad::pop(_t72);
                    _d_cur_nnz -= _r_d16;
                }
            }
            clad::pop(_t71);
        }
        {
            double _r_d15 = _d_normr;
            double _r37 = _r_d15 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t70), 1.).pushforward;
            _d_rtrans += _r37;
            double _r38 = clad::pop(_EERepl_normr3);
            _delta_normr += clad::getErrorVal(_r_d15, _r38, "normr");
            _d_normr -= _r_d15;
        }
        if (clad::pop(_t45)) {
            {
                for (; clad::back(_t46); clad::back(_t46)--) {
                    int _t48 = clad::pop(_t47);
                    double _r_d10 = _d_p[_t48];
                    int _t50 = clad::pop(_t49);
                    _d_r[_t50] += _r_d10;
                    double _r19 = _r_d10 * clad::pop(_t51);
                    _d_beta += _r19;
                    double _r20 = clad::pop(_t52) * _r_d10;
                    int _t54 = clad::pop(_t53);
                    _d_r[_t54] += _r20;
                    double _r21 = clad::pop(_EERepl_p2);
                    int _r22 = clad::pop(_t47);
                    _delta_p[_r22] += clad::getErrorVal(_r_d10, _r21, "p");
                    _final_error += _delta_p[_r22];
                    _d_p[_t48] -= _r_d10;
                    _d_p[_t48];
                }
                clad::pop(_t46);
            }
            {
                double _r_d9 = _d_beta;
                double _r18 = clad::pop(_EERepl_beta2);
                _delta_beta += clad::getErrorVal(_r_d9, _r18, "beta");
                _d_beta -= _r_d9;
            }
        } else {
            {
                for (; clad::back(_t61); clad::back(_t61)--) {
                    int _t63 = clad::pop(_t62);
                    double _r_d14 = _d_p[_t63];
                    int _t65 = clad::pop(_t64);
                    _d_r[_t65] += _r_d14;
                    double _r33 = _r_d14 * clad::pop(_t66);
                    _d_beta += _r33;
                    double _r34 = clad::pop(_t67) * _r_d14;
                    int _t69 = clad::pop(_t68);
                    _d_p[_t69] += _r34;
                    double _r35 = clad::pop(_EERepl_p3);
                    int _r36 = clad::pop(_t62);
                    _delta_p[_r36] += clad::getErrorVal(_r_d14, _r35, "p");
                    _final_error += _delta_p[_r36];
                    _d_p[_t63] -= _r_d14;
                    _d_p[_t63];
                }
                clad::pop(_t61);
            }
            {
                double _r_d13 = _d_beta;
                double _r29 = clad::pop(_t59);
                double _r30 = _r_d13 / _r29;
                _d_rtrans += _r30;
                double _r31 = _r_d13 * -clad::pop(_t60) / (_r29 * _r29);
                _d_oldrtrans += _r31;
                double _r32 = clad::pop(_EERepl_beta3);
                _delta_beta += clad::getErrorVal(_r_d13, _r32, "beta");
                _d_beta -= _r_d13;
            }
            {
                double _r_d12 = _d_rtrans;
                double *_r25 = clad::pop(_t56);
                double *_r27 = clad::pop(_t57);
                int _grad3 = 0;
                double _t58 = 0;
                ddot_pullback(clad::pop(_t55), _r25, _r27, _r_d12, &_grad3, _d_r, _d_r, _t58);
                int _r23 = _grad3;
                _d_nrow += _r23;
                clad::array<double> _r24(_d_r);
                clad::array<double> _r26(_d_r);
                double _r28 = clad::pop(_EERepl_rtrans2);
                _delta_rtrans += _t58;
                _d_rtrans -= _r_d12;
            }
            {
                double _r_d11 = _d_oldrtrans;
                _d_rtrans += _r_d11;
                _d_oldrtrans -= _r_d11;
            }
        }
    }
    {
        double _r_d8 = _d_normr;
        double _r17 = _r_d8 * clad::custom_derivatives::sqrt_pushforward(_t42, 1.).pushforward;
        _d_rtrans += _r17;
        _delta_normr += clad::getErrorVal(_r_d8, _EERepl_normr2, "normr");
        _d_normr -= _r_d8;
    }
    {
        double _r_d7 = _d_rtrans;
        int _grad0 = 0;
        double _t41 = 0;
        ddot_pullback(_t38, _t39, _t40, _r_d7, &_grad0, _d_r, _d_r, _t41);
        int _r14 = _grad0;
        _d_nrow += _r14;
        clad::array<double> _r15(_d_r);
        clad::array<double> _r16(_d_r);
        _delta_rtrans += _t41;
        _d_rtrans -= _r_d7;
    }
    for (; _t29; _t29--) {
        int _t31 = clad::pop(_t30);
        double _r_d6 = _d_r[_t31];
        int _t33 = clad::pop(_t32);
        _d_b[_t33] += _r_d6;
        double _r10 = _r_d6 * clad::pop(_t34);
        _d_beta += _r10;
        double _r11 = clad::pop(_t35) * _r_d6;
        int _t37 = clad::pop(_t36);
        _d_Ap[_t37] += _r11;
        double _r12 = clad::pop(_EERepl_r1);
        int _r13 = clad::pop(_t30);
        _delta_r[_r13] += clad::getErrorVal(_r_d6, _r12, "r");
        _final_error += _delta_r[_r13];
        _d_r[_t31] -= _r_d6;
        _d_r[_t31];
    }
    {
        int _t26 = clad::pop(_t25);
        int _t24 = clad::pop(_t23);
        double _r_d5 = _d_beta;
        _delta_beta += clad::getErrorVal(_r_d5, _EERepl_beta1, "beta");
        _d_beta -= _r_d5;
    }
    for (; _t9; _t9--) {
        {
            for (; clad::back(_t14); clad::back(_t14)--) {
                int _t16 = clad::pop(_t15);
                double _r_d4 = _d_Ap[_t16];
                _d_Ap[_t16] += _r_d4;
                double _r6 = _r_d4 * clad::pop(_t17);
                int _t19 = clad::pop(_t18);
                int _t21 = clad::pop(_t20);
                double _r7 = clad::pop(_t22) * _r_d4;
                int _t28 = clad::pop(_t27);
                _d_p[_t28] += _r7;
                double _r8 = clad::pop(_EERepl_Ap2);
                int _r9 = clad::pop(_t15);
                _delta_Ap[_r9] += clad::getErrorVal(_r_d4, _r8, "Ap");
                _final_error += _delta_Ap[_r9];
                _d_Ap[_t16] -= _r_d4;
                _d_Ap[_t16];
            }
            clad::pop(_t14);
        }
        {
            int _t13 = clad::pop(_t12);
            double _r_d3 = _d_Ap[_t13];
            double _r4 = clad::pop(_EERepl_Ap1);
            int _r5 = clad::pop(_t12);
            _delta_Ap[_r5] += clad::getErrorVal(_r_d3, _r4, "Ap");
            _final_error += _delta_Ap[_r5];
            _d_Ap[_t13] -= _r_d3;
            _d_Ap[_t13];
        }
        {
            int _r_d2 = _d_cur_nnz;
            int _t11 = clad::pop(_t10);
            _d_cur_nnz -= _r_d2;
        }
    }
    for (; _t0; _t0--) {
        int _t2 = clad::pop(_t1);
        double _r_d1 = _d_p[_t2];
        int _t4 = clad::pop(_t3);
        _d_x[_t4] += _r_d1;
        double _r0 = _r_d1 * clad::pop(_t5);
        _d_beta += _r0;
        double _r1 = clad::pop(_t6) * _r_d1;
        int _t8 = clad::pop(_t7);
        _d_x[_t8] += _r1;
        double _r2 = clad::pop(_EERepl_p1);
        int _r3 = clad::pop(_t1);
        _delta_p[_r3] += clad::getErrorVal(_r_d1, _r2, "p");
        _final_error += _delta_p[_r3];
        _d_p[_t2] -= _r_d1;
        _d_p[_t2];
    }
    _delta_beta += clad::getErrorVal(_d_beta, _EERepl_beta0, "beta");
    _delta_oldrtrans += clad::getErrorVal(_d_oldrtrans, _EERepl_oldrtrans0, "oldrtrans");
    _delta_rtrans += clad::getErrorVal(_d_rtrans, _EERepl_rtrans0, "rtrans");
    {
        double _r_d0 = _d_normr;
        _delta_normr += clad::getErrorVal(_r_d0, _EERepl_normr1, "normr");
        _d_normr -= _r_d0;
    }
    _delta_tolerance += clad::getErrorVal(_d_tolerance, _EERepl_tolerance0, "tolerance");
    _delta_normr += clad::getErrorVal(_d_normr, _EERepl_normr0, "normr");
    clad::array<double> _delta_b(_d_b.size());
    int i = 0;
    for (; i < _d_b.size(); i++) {
        double _t124 = clad::getErrorVal(_d_b[i], b[i], "b");
        _delta_b[i] += _t124;
        _final_error += _t124;
    }
    i = 0;
    for (; i < _d_x.size(); i++) {
        double _t125 = clad::getErrorVal(_d_x[i], _EERepl_x0[i], "x");
        _delta_x[i] += _t125;
        _final_error += _t125;
    }
    clad::array<double> _delta_xexact(_d_xexact.size());
    i = 0;
    for (; i < _d_xexact.size(); i++) {
        double _t126 = clad::getErrorVal(_d_xexact[i], xexact[i], "xexact");
        _delta_xexact[i] += _t126;
        _final_error += _t126;
    }
    i = 0;
    for (; i < _d_r.size(); i++) {
        double _t127 = clad::getErrorVal(_d_r[i], _EERepl_r0[i], "r");
        _delta_r[i] += _t127;
        _final_error += _t127;
    }
    i = 0;
    for (; i < _d_p.size(); i++) {
        double _t128 = clad::getErrorVal(_d_p[i], _EERepl_p0[i], "p");
        _delta_p[i] += _t128;
        _final_error += _t128;
    }
    i = 0;
    for (; i < _d_Ap.size(); i++) {
        double _t129 = clad::getErrorVal(_d_Ap[i], _EERepl_Ap0[i], "Ap");
        _delta_Ap[i] += _t129;
        _final_error += _t129;
    }
    _final_error += _delta_residual + _delta_rtrans + _delta_oldrtrans + _delta_tolerance + _delta_diff + _delta_beta + _delta_normr;
}
