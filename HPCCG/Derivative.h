void ddot_pullback(int n, double *x, double *y, double _d_y0, clad::array_ref<int> _d_n, clad::array_ref<double> _d_x, clad::array_ref<double> _d_y, double &_final_error) {
    double _delta_local_result = 0;
    double _EERepl_local_result0;
    bool _cond0;
    unsigned long _t0;
    clad::tape<double> _t1 = {};
    clad::tape<int> _t2 = {};
    clad::tape<double> _t4 = {};
    clad::tape<int> _t5 = {};
    clad::tape<double> _EERepl_local_result1 = {};
    unsigned long _t7;
    clad::tape<double> _t8 = {};
    clad::tape<int> _t9 = {};
    clad::tape<double> _t11 = {};
    clad::tape<int> _t12 = {};
    clad::tape<double> _EERepl_local_result2 = {};
    double _d_local_result = 0;
    double local_result = 0.;
    _EERepl_local_result0 = local_result;
    _cond0 = y == x;
    if (_cond0) {
        _t0 = 0;
        {
            int _d_i = 0;
            for (int i = 0; i < n; i++) {
                _t0++;
                local_result += clad::push(_t4, x[clad::push(_t2, i)]) * clad::push(_t1, x[clad::push(_t5, i)]);
                clad::push(_EERepl_local_result1, local_result);
            }
        }
    } else {
        _t7 = 0;
        {
            int _d_i = 0;
            for (int i = 0; i < n; i++) {
                _t7++;
                local_result += clad::push(_t11, x[clad::push(_t9, i)]) * clad::push(_t8, y[clad::push(_t12, i)]);
                clad::push(_EERepl_local_result2, local_result);
            }
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
            _delta_local_result += std::abs(_r_d0 * _r2 * 1.1920928955078125E-7);
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
            _delta_local_result += std::abs(_r_d1 * _r5 * 1.1920928955078125E-7);
            _d_local_result -= _r_d1;
        }
    _delta_local_result += std::abs(_d_local_result * _EERepl_local_result0 * 1.1920928955078125E-7);
    clad::array<double> _delta_x(_d_x.size());
    int i = 0;
    for (; i < _d_x.size(); i++) {
        double _t14 = std::abs(_d_x[i] * x[i] * 1.1920928955078125E-7);
        _delta_x[i] += _t14;
        _final_error += _t14;
    }
    clad::array<double> _delta_y(_d_y.size());
    i = 0;
    for (; i < _d_y.size(); i++) {
        double _t15 = std::abs(_d_y[i] * y[i] * 1.1920928955078125E-7);
        _delta_y[i] += _t15;
        _final_error += _t15;
    }
    _final_error += _delta_local_result;
}

void HPCCG_residual_grad(double *b, double *x, double *xexact, double *r, double *p, double *Ap, clad::array_ref<double> _d_b, clad::array_ref<double> _d_x, clad::array_ref<double> _d_xexact, clad::array_ref<double> _d_r, clad::array_ref<double> _d_p, clad::array_ref<double> _d_Ap, double &_final_error) {
    double _delta_normr = 0;
    double _EERepl_normr0;
    double _delta_tolerance = 0;
    double _EERepl_tolerance0;
    double _delta_sum = 0;
    double _EERepl_sum0;
    double _EERepl_normr1;
    double _delta_rtrans = 0;
    double _EERepl_rtrans0;
    double _delta_oldrtrans = 0;
    double _EERepl_oldrtrans0;
    double _delta_alpha = 0;
    double _EERepl_alpha0;
    double _delta_alpha0 = 0;
    double _EERepl_alpha1;
    unsigned long _t0;
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
    clad::tape<double> _EERepl_sum1 = {};
    clad::tape<int> _t10 = {};
    clad::tape<unsigned long> _t12 = {};
    clad::tape<double> _t13 = {};
    clad::tape<int> _t14 = {};
    clad::tape<int> _t16 = {};
    clad::tape<double> _t18 = {};
    clad::tape<int> _t19 = {};
    clad::tape<int> _t21 = {};
    clad::tape<int> _t23 = {};
    clad::tape<double> _EERepl_sum2 = {};
    clad::tape<int> _t25 = {};
    double _delta_beta = 0;
    double _EERepl_beta0;
    double _EERepl_beta1;
    unsigned long _t27;
    clad::tape<int> _t28 = {};
    clad::array<double> _delta_r(_d_r.size());
    clad::array<double> _EERepl_r0(_d_r.size());
    for (int i0 = 0; i0 < _d_r.size(); i0++) {
        _EERepl_r0[i0] = r[i0];
    }
    clad::tape<int> _t30 = {};
    clad::tape<double> _t32 = {};
    clad::tape<double> _t33 = {};
    clad::tape<int> _t34 = {};
    clad::tape<double> _EERepl_r1 = {};
    int _t36;
    double *_t37;
    double *_t38;
    double _EERepl_rtrans1;
    double _t40;
    double _EERepl_normr2;
    unsigned long _t41;
    clad::tape<bool> _t43 = {};
    clad::tape<double> _EERepl_beta2 = {};
    clad::tape<unsigned long> _t44 = {};
    clad::tape<int> _t45 = {};
    clad::tape<int> _t47 = {};
    clad::tape<double> _t49 = {};
    clad::tape<double> _t50 = {};
    clad::tape<int> _t51 = {};
    clad::tape<double> _EERepl_p2 = {};
    clad::tape<int> _t53 = {};
    clad::tape<double *> _t54 = {};
    clad::tape<double *> _t55 = {};
    clad::tape<double> _EERepl_rtrans2 = {};
    clad::tape<double> _t57 = {};
    clad::tape<double> _t58 = {};
    clad::tape<double> _EERepl_beta3 = {};
    clad::tape<unsigned long> _t59 = {};
    clad::tape<int> _t60 = {};
    clad::tape<int> _t62 = {};
    clad::tape<double> _t64 = {};
    clad::tape<double> _t65 = {};
    clad::tape<int> _t66 = {};
    clad::tape<double> _EERepl_p3 = {};
    clad::tape<double> _t68 = {};
    clad::tape<double> _EERepl_normr3 = {};
    clad::tape<unsigned long> _t69 = {};
    clad::tape<double> _EERepl_sum3 = {};
    clad::tape<int> _t70 = {};
    clad::tape<unsigned long> _t72 = {};
    clad::tape<double> _t73 = {};
    clad::tape<int> _t74 = {};
    clad::tape<int> _t76 = {};
    clad::tape<double> _t78 = {};
    clad::tape<int> _t79 = {};
    clad::tape<int> _t81 = {};
    clad::tape<int> _t83 = {};
    clad::tape<double> _EERepl_sum4 = {};
    clad::tape<int> _t85 = {};
    clad::tape<int> _t87 = {};
    clad::tape<double *> _t88 = {};
    clad::tape<double *> _t89 = {};
    clad::tape<double> _EERepl_beta4 = {};
    clad::tape<double> _t91 = {};
    clad::tape<double> _t92 = {};
    clad::tape<double> _EERepl_beta5 = {};
    clad::tape<unsigned long> _t93 = {};
    clad::tape<int> _t94 = {};
    clad::array<double> _delta_x(_d_x.size());
    clad::array<double> _EERepl_x0(_d_x.size());
    for (int i0 = 0; i0 < _d_x.size(); i0++) {
        _EERepl_x0[i0] = x[i0];
    }
    clad::tape<int> _t96 = {};
    clad::tape<double> _t98 = {};
    clad::tape<double> _t99 = {};
    clad::tape<int> _t100 = {};
    clad::tape<double> _EERepl_x1 = {};
    clad::tape<double> _EERepl_beta6 = {};
    clad::tape<unsigned long> _t102 = {};
    clad::tape<int> _t103 = {};
    clad::tape<int> _t105 = {};
    clad::tape<double> _t107 = {};
    clad::tape<double> _t108 = {};
    clad::tape<int> _t109 = {};
    clad::tape<double> _EERepl_r2 = {};
    double _delta_residual = 0;
    double _EERepl_residual0;
    double _delta_diff = 0;
    double _EERepl_diff0;
    unsigned long _t111;
    clad::tape<int> _t112 = {};
    clad::tape<int> _t114 = {};
    clad::tape<double> _t116 = {};
    clad::tape<double> _EERepl_diff1 = {};
    clad::tape<bool> _t119 = {};
    int _d_niters = 0;
    int niters = 0;
    double _d_normr = 0;
    double normr = 0.;
    _EERepl_normr0 = normr;
    int _d_max_iter = 0;
    int max_iter = 100;
    double _d_tolerance = 0;
    double tolerance = 0.;
    _EERepl_tolerance0 = tolerance;
    double _d_sum = 0;
    double sum;
    _EERepl_sum0 = sum;
    int _d_cur_nnz = 0;
    int cur_nnz;
    int _d_nrow = 0;
    int nrow = A.local_nrow;
    int _d_ncol = 0;
    int ncol = A.local_ncol;
    normr = 0.;
    _EERepl_normr1 = normr;
    double _d_rtrans = 0;
    double rtrans = 0.;
    _EERepl_rtrans0 = rtrans;
    double _d_oldrtrans = 0;
    double oldrtrans = 0.;
    _EERepl_oldrtrans0 = oldrtrans;
    double _d_alpha = 0, _d_beta = 0;
    double alpha = 1., beta = 0.;
    _EERepl_beta0 = beta;
    _EERepl_alpha1 = alpha;
    _EERepl_alpha0 = alpha;
    _t0 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < nrow; i++) {
            _t0++;
            p[clad::push(_t1, i)] = x[clad::push(_t3, i)] + clad::push(_t6, beta) * clad::push(_t5, x[clad::push(_t7, i)]);
            clad::push(_EERepl_p1, p[clad::push(_t1, i)]);
        }
    }
    _t9 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < nrow; i++) {
            _t9++;
            sum = 0.;
            clad::push(_EERepl_sum1, sum);
            cur_nnz = A.nnz_in_row[clad::push(_t10, i)];
            clad::push(_t12, 0UL);
            {
                int _d_j = 0;
                for (int j = 0; j < cur_nnz; j++) {
                    clad::back(_t12)++;
                    sum += clad::push(_t18, A.ptr_to_vals_in_row[clad::push(_t14, i)][clad::push(_t16, j)]) * clad::push(_t13, p[clad::push(_t23, A.ptr_to_inds_in_row[clad::push(_t19, i)][clad::push(_t21, j)])]);
                    clad::push(_EERepl_sum2, sum);
                }
            }
            Ap[clad::push(_t25, i)] = sum;
        }
    }
    beta = -1.;
    _EERepl_beta1 = beta;
    _t27 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < nrow; i++) {
            _t27++;
            r[clad::push(_t28, i)] = b[clad::push(_t30, i)] + clad::push(_t33, beta) * clad::push(_t32, Ap[clad::push(_t34, i)]);
            clad::push(_EERepl_r1, r[clad::push(_t28, i)]);
        }
    }
    _t36 = nrow;
    _t37 = r;
    _t38 = r;
    rtrans = ddot(_t36, r, r);
    _EERepl_rtrans1 = rtrans;
    _t40 = rtrans;
    normr = sqrt(_t40);
    _EERepl_normr2 = normr;
    _t41 = 0;
    int _d_k = 0;
    {
        _d_k = 0;
        for (int k = 1; k < max_iter && normr > tolerance; k++) {
            _t41++;
            bool _t42 = k == 1;
            {
                if (_t42) {
                    beta = 0.;
                    clad::push(_EERepl_beta2, beta);
                    clad::push(_t44, 0UL);
                    {
                        int _d_i = 0;
                        for (int i = 0; i < nrow; i++) {
                            clad::back(_t44)++;
                            p[clad::push(_t45, i)] = r[clad::push(_t47, i)] + clad::push(_t50, beta) * clad::push(_t49, r[clad::push(_t51, i)]);
                            clad::push(_EERepl_p2, p[clad::push(_t45, i)]);
                        }
                    }
                } else {
                    oldrtrans = rtrans;
                    clad::push(_t54, r);
                    clad::push(_t55, r);
                    rtrans = ddot(clad::push(_t53, nrow), r, r);
                    clad::push(_EERepl_rtrans2, rtrans);
                    beta = clad::push(_t58, rtrans) / clad::push(_t57, oldrtrans);
                    clad::push(_EERepl_beta3, beta);
                    clad::push(_t59, 0UL);
                    {
                        int _d_i = 0;
                        for (int i = 0; i < nrow; i++) {
                            clad::back(_t59)++;
                            p[clad::push(_t60, i)] = r[clad::push(_t62, i)] + clad::push(_t65, beta) * clad::push(_t64, p[clad::push(_t66, i)]);
                            clad::push(_EERepl_p3, p[clad::push(_t60, i)]);
                        }
                    }
                }
                clad::push(_t43, _t42);
            }
            normr = sqrt(clad::push(_t68, rtrans));
            clad::push(_EERepl_normr3, normr);
            clad::push(_t69, 0UL);
            {
                int _d_i = 0;
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t69)++;
                    sum = 0.;
                    clad::push(_EERepl_sum3, sum);
                    cur_nnz = A.nnz_in_row[clad::push(_t70, i)];
                    clad::push(_t72, 0UL);
                    {
                        int _d_j = 0;
                        for (int j = 0; j < cur_nnz; j++) {
                            clad::back(_t72)++;
                            sum += clad::push(_t78, A.ptr_to_vals_in_row[clad::push(_t74, i)][clad::push(_t76, j)]) * clad::push(_t73, p[clad::push(_t83, A.ptr_to_inds_in_row[clad::push(_t79, i)][clad::push(_t81, j)])]);
                            clad::push(_EERepl_sum4, sum);
                        }
                    }
                    Ap[clad::push(_t85, i)] = sum;
                }
            }
            clad::push(_t88, p);
            clad::push(_t89, Ap);
            beta = ddot(clad::push(_t87, nrow), p, Ap);
            clad::push(_EERepl_beta4, beta);
            beta = clad::push(_t92, rtrans) / clad::push(_t91, beta);
            clad::push(_EERepl_beta5, beta);
            clad::push(_t93, 0UL);
            {
                int _d_i = 0;
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t93)++;
                    x[clad::push(_t94, i)] = x[clad::push(_t96, i)] + clad::push(_t99, beta) * clad::push(_t98, p[clad::push(_t100, i)]);
                    clad::push(_EERepl_x1, x[clad::push(_t94, i)]);
                }
            }
            beta = -beta;
            clad::push(_EERepl_beta6, beta);
            clad::push(_t102, 0UL);
            {
                int _d_i = 0;
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t102)++;
                    r[clad::push(_t103, i)] = r[clad::push(_t105, i)] + clad::push(_t108, beta) * clad::push(_t107, Ap[clad::push(_t109, i)]);
                    clad::push(_EERepl_r2, r[clad::push(_t103, i)]);
                }
            }
            niters = k;
        }
    }
    double _d_residual = 0;
    double residual = 0.;
    _EERepl_residual0 = residual;
    double _d_diff = 0;
    double diff;
    _EERepl_diff0 = diff;
    _t111 = 0;
    {
        int _d_i = 0;
        for (int i = 0; i < nrow; i++) {
            _t111++;
            diff = fabs(clad::push(_t116, x[clad::push(_t112, i)] - xexact[clad::push(_t114, i)]));
            clad::push(_EERepl_diff1, diff);
            bool _t118 = diff > residual;
            {
                if (_t118)
                    residual = diff;
                clad::push(_t119, _t118);
            }
        }
    }
    double HPCCG_residual_return = residual;
    goto _label0;
  _label0:
    _d_residual += 1;
    for (; _t111; _t111--) {
        if (clad::pop(_t119)) {
            double _r_d28 = _d_residual;
            _d_diff += _r_d28;
            _d_residual -= _r_d28;
        }
        {
            double _r_d27 = _d_diff;
            double _changed_ = clad::pop(_t116);
            double _r60 =
                _r_d27 * numerical_diff::forward_central_difference(
                             fabs, _changed_, 0, 0, _changed_);
            int _t113 = clad::pop(_t112);
            _d_x[_t113] += _r60;
            int _t115 = clad::pop(_t114);
            _d_xexact[_t115] += -_r60;
            double _r61 = clad::pop(_EERepl_diff1);
            _delta_diff += std::abs(_r_d27 * _r61 * 1.1920928955078125E-7);
            _d_diff -= _r_d27;
        }
    }
    _delta_residual += std::abs(_d_residual * _EERepl_residual0 * 1.1920928955078125E-7);
    for (; _t41; _t41--) {
        {
            int _r_d26 = _d_niters;
            _d_k += _r_d26;
            _d_niters -= _r_d26;
        }
        {
            for (; clad::back(_t102); clad::back(_t102)--) {
                int _t104 = clad::pop(_t103);
                double _r_d25 = _d_r[_t104];
                int _t106 = clad::pop(_t105);
                _d_r[_t106] += _r_d25;
                double _r56 = _r_d25 * clad::pop(_t107);
                _d_beta += _r56;
                double _r57 = clad::pop(_t108) * _r_d25;
                int _t110 = clad::pop(_t109);
                _d_Ap[_t110] += _r57;
                double _r58 = clad::pop(_EERepl_r2);
                int _r59 = clad::pop(_t103);
                _delta_r[_r59] += std::abs(_r_d25 * _r58 * 1.1920928955078125E-7);
                _final_error += _delta_r[_r59];
                _d_r[_t104] -= _r_d25;
                _d_r[_t104];
            }
            clad::pop(_t102);
        }
        {
            double _r_d24 = _d_beta;
            _d_beta += -_r_d24;
            double _r55 = clad::pop(_EERepl_beta6);
            _delta_beta += std::abs(_r_d24 * _r55 * 1.1920928955078125E-7);
            _d_beta -= _r_d24;
        }
        {
            for (; clad::back(_t93); clad::back(_t93)--) {
                int _t95 = clad::pop(_t94);
                double _r_d23 = _d_x[_t95];
                int _t97 = clad::pop(_t96);
                _d_x[_t97] += _r_d23;
                double _r51 = _r_d23 * clad::pop(_t98);
                _d_beta += _r51;
                double _r52 = clad::pop(_t99) * _r_d23;
                int _t101 = clad::pop(_t100);
                _d_p[_t101] += _r52;
                double _r53 = clad::pop(_EERepl_x1);
                int _r54 = clad::pop(_t94);
                _delta_x[_r54] += std::abs(_r_d23 * _r53 * 1.1920928955078125E-7);
                _final_error += _delta_x[_r54];
                _d_x[_t95] -= _r_d23;
                _d_x[_t95];
            }
            clad::pop(_t93);
        }
        {
            double _r_d22 = _d_beta;
            double _r47 = clad::pop(_t91);
            double _r48 = _r_d22 / _r47;
            _d_rtrans += _r48;
            double _r49 = _r_d22 * -clad::pop(_t92) / (_r47 * _r47);
            _d_beta += _r49;
            double _r50 = clad::pop(_EERepl_beta5);
            _delta_beta += std::abs(_r_d22 * _r50 * 1.1920928955078125E-7);
            _d_beta -= _r_d22;
        }
        {
            double _r_d21 = _d_beta;
            double *_r43 = clad::pop(_t88);
            double *_r45 = clad::pop(_t89);
            int _grad6 = 0;
            double _t90 = 0;
            ddot_pullback(clad::pop(_t87), _r43, _r45, _r_d21, &_grad6, _d_p, _d_Ap, _t90);
            int _r41 = _grad6;
            _d_nrow += _r41;
            clad::array<double> _r42(_d_p);
            clad::array<double> _r44(_d_Ap);
            double _r46 = clad::pop(_EERepl_beta4);
            _delta_beta += _t90;
            _d_beta -= _r_d21;
        }
        {
            for (; clad::back(_t69); clad::back(_t69)--) {
                {
                    int _t86 = clad::pop(_t85);
                    int _t82 = clad::pop(_t81);
                    int _t80 = clad::pop(_t79);
                    double _r_d20 = _d_Ap[_t86];
                    _d_sum += _r_d20;
                    _d_Ap[_t86] -= _r_d20;
                    _d_Ap[_t86];
                }
                {
                    for (; clad::back(_t72); clad::back(_t72)--) {
                        double _r_d19 = _d_sum;
                        _d_sum += _r_d19;
                        double _r38 = _r_d19 * clad::pop(_t73);
                        int _t75 = clad::pop(_t74);
                        int _t77 = clad::pop(_t76);
                        double _r39 = clad::pop(_t78) * _r_d19;
                        int _t84 = clad::pop(_t83);
                        _d_p[_t84] += _r39;
                        double _r40 = clad::pop(_EERepl_sum4);
                        _delta_sum += std::abs(_r_d19 * _r40 * 1.1920928955078125E-7);
                        _d_sum -= _r_d19;
                    }
                    clad::pop(_t72);
                }
                {
                    int _r_d18 = _d_cur_nnz;
                    int _t71 = clad::pop(_t70);
                    _d_cur_nnz -= _r_d18;
                }
                {
                    double _r_d17 = _d_sum;
                    double _r37 = clad::pop(_EERepl_sum3);
                    _delta_sum += std::abs(_r_d17 * _r37 * 1.1920928955078125E-7);
                    _d_sum -= _r_d17;
                }
            }
            clad::pop(_t69);
        }
        {
            double _r_d16 = _d_normr;
            double _r35 = _r_d16 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t68), 1.).pushforward;
            _d_rtrans += _r35;
            double _r36 = clad::pop(_EERepl_normr3);
            _delta_normr += std::abs(_r_d16 * _r36 * 1.1920928955078125E-7);
            _d_normr -= _r_d16;
        }
        if (clad::pop(_t43)) {
            {
                for (; clad::back(_t44); clad::back(_t44)--) {
                    int _t46 = clad::pop(_t45);
                    double _r_d11 = _d_p[_t46];
                    int _t48 = clad::pop(_t47);
                    _d_r[_t48] += _r_d11;
                    double _r17 = _r_d11 * clad::pop(_t49);
                    _d_beta += _r17;
                    double _r18 = clad::pop(_t50) * _r_d11;
                    int _t52 = clad::pop(_t51);
                    _d_r[_t52] += _r18;
                    double _r19 = clad::pop(_EERepl_p2);
                    int _r20 = clad::pop(_t45);
                    _delta_p[_r20] += std::abs(_r_d11 * _r19 * 1.1920928955078125E-7);
                    _final_error += _delta_p[_r20];
                    _d_p[_t46] -= _r_d11;
                    _d_p[_t46];
                }
                clad::pop(_t44);
            }
            {
                double _r_d10 = _d_beta;
                double _r16 = clad::pop(_EERepl_beta2);
                _delta_beta += std::abs(_r_d10 * _r16 * 1.1920928955078125E-7);
                _d_beta -= _r_d10;
            }
        } else {
            {
                for (; clad::back(_t59); clad::back(_t59)--) {
                    int _t61 = clad::pop(_t60);
                    double _r_d15 = _d_p[_t61];
                    int _t63 = clad::pop(_t62);
                    _d_r[_t63] += _r_d15;
                    double _r31 = _r_d15 * clad::pop(_t64);
                    _d_beta += _r31;
                    double _r32 = clad::pop(_t65) * _r_d15;
                    int _t67 = clad::pop(_t66);
                    _d_p[_t67] += _r32;
                    double _r33 = clad::pop(_EERepl_p3);
                    int _r34 = clad::pop(_t60);
                    _delta_p[_r34] += std::abs(_r_d15 * _r33 * 1.1920928955078125E-7);
                    _final_error += _delta_p[_r34];
                    _d_p[_t61] -= _r_d15;
                    _d_p[_t61];
                }
                clad::pop(_t59);
            }
            {
                double _r_d14 = _d_beta;
                double _r27 = clad::pop(_t57);
                double _r28 = _r_d14 / _r27;
                _d_rtrans += _r28;
                double _r29 = _r_d14 * -clad::pop(_t58) / (_r27 * _r27);
                _d_oldrtrans += _r29;
                double _r30 = clad::pop(_EERepl_beta3);
                _delta_beta += std::abs(_r_d14 * _r30 * 1.1920928955078125E-7);
                _d_beta -= _r_d14;
            }
            {
                double _r_d13 = _d_rtrans;
                double *_r23 = clad::pop(_t54);
                double *_r25 = clad::pop(_t55);
                int _grad3 = 0;
                double _t56 = 0;
                ddot_pullback(clad::pop(_t53), _r23, _r25, _r_d13, &_grad3, _d_r, _d_r, _t56);
                int _r21 = _grad3;
                _d_nrow += _r21;
                clad::array<double> _r22(_d_r);
                clad::array<double> _r24(_d_r);
                double _r26 = clad::pop(_EERepl_rtrans2);
                _delta_rtrans += _t56;
                _d_rtrans -= _r_d13;
            }
            {
                double _r_d12 = _d_oldrtrans;
                _d_rtrans += _r_d12;
                _d_oldrtrans -= _r_d12;
            }
        }
    }
    {
        double _r_d9 = _d_normr;
        double _r15 = _r_d9 * clad::custom_derivatives::sqrt_pushforward(_t40, 1.).pushforward;
        _d_rtrans += _r15;
        _delta_normr += std::abs(_r_d9 * _EERepl_normr2 * 1.1920928955078125E-7);
        _d_normr -= _r_d9;
    }
    {
        double _r_d8 = _d_rtrans;
        int _grad0 = 0;
        double _t39 = 0;
        ddot_pullback(_t36, _t37, _t38, _r_d8, &_grad0, _d_r, _d_r, _t39);
        int _r12 = _grad0;
        _d_nrow += _r12;
        clad::array<double> _r13(_d_r);
        clad::array<double> _r14(_d_r);
        _delta_rtrans += _t39;
        _d_rtrans -= _r_d8;
    }
    for (; _t27; _t27--) {
        int _t29 = clad::pop(_t28);
        double _r_d7 = _d_r[_t29];
        int _t31 = clad::pop(_t30);
        _d_b[_t31] += _r_d7;
        double _r8 = _r_d7 * clad::pop(_t32);
        _d_beta += _r8;
        double _r9 = clad::pop(_t33) * _r_d7;
        int _t35 = clad::pop(_t34);
        _d_Ap[_t35] += _r9;
        double _r10 = clad::pop(_EERepl_r1);
        int _r11 = clad::pop(_t28);
        _delta_r[_r11] += std::abs(_r_d7 * _r10 * 1.1920928955078125E-7);
        _final_error += _delta_r[_r11];
        _d_r[_t29] -= _r_d7;
        _d_r[_t29];
    }
    {
        double _r_d6 = _d_beta;
        _delta_beta += std::abs(_r_d6 * _EERepl_beta1 * 1.1920928955078125E-7);
        _d_beta -= _r_d6;
    }
    for (; _t9; _t9--) {
        {
            int _t26 = clad::pop(_t25);
            int _t22 = clad::pop(_t21);
            int _t20 = clad::pop(_t19);
            double _r_d5 = _d_Ap[_t26];
            _d_sum += _r_d5;
            _d_Ap[_t26] -= _r_d5;
            _d_Ap[_t26];
        }
        {
            for (; clad::back(_t12); clad::back(_t12)--) {
                double _r_d4 = _d_sum;
                _d_sum += _r_d4;
                double _r5 = _r_d4 * clad::pop(_t13);
                int _t15 = clad::pop(_t14);
                int _t17 = clad::pop(_t16);
                double _r6 = clad::pop(_t18) * _r_d4;
                int _t24 = clad::pop(_t23);
                _d_p[_t24] += _r6;
                double _r7 = clad::pop(_EERepl_sum2);
                _delta_sum += std::abs(_r_d4 * _r7 * 1.1920928955078125E-7);
                _d_sum -= _r_d4;
            }
            clad::pop(_t12);
        }
        {
            int _r_d3 = _d_cur_nnz;
            int _t11 = clad::pop(_t10);
            _d_cur_nnz -= _r_d3;
        }
        {
            double _r_d2 = _d_sum;
            double _r4 = clad::pop(_EERepl_sum1);
            _delta_sum += std::abs(_r_d2 * _r4 * 1.1920928955078125E-7);
            _d_sum -= _r_d2;
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
        _delta_p[_r3] += std::abs(_r_d1 * _r2 * 1.1920928955078125E-7);
        _final_error += _delta_p[_r3];
        _d_p[_t2] -= _r_d1;
        _d_p[_t2];
    }
    {
        _delta_alpha += std::abs(_d_alpha * _EERepl_alpha0 * 1.1920928955078125E-7);
        _delta_alpha0 += std::abs(_d_alpha * _EERepl_alpha1 * 1.1920928955078125E-7);
    }
    _delta_oldrtrans += std::abs(_d_oldrtrans * _EERepl_oldrtrans0 * 1.1920928955078125E-7);
    _delta_rtrans += std::abs(_d_rtrans * _EERepl_rtrans0 * 1.1920928955078125E-7);
    {
        double _r_d0 = _d_normr;
        _delta_normr += std::abs(_r_d0 * _EERepl_normr1 * 1.1920928955078125E-7);
        _d_normr -= _r_d0;
    }
    _delta_tolerance += std::abs(_d_tolerance * _EERepl_tolerance0 * 1.1920928955078125E-7);
    _delta_normr += std::abs(_d_normr * _EERepl_normr0 * 1.1920928955078125E-7);
    clad::array<double> _delta_b(_d_b.size());
    int i = 0;
    for (; i < _d_b.size(); i++) {
        double _t120 = std::abs(_d_b[i] * b[i] * 1.1920928955078125E-7);
        _delta_b[i] += _t120;
        _final_error += _t120;
    }
    i = 0;
    for (; i < _d_x.size(); i++) {
        double _t121 = std::abs(_d_x[i] * _EERepl_x0[i] * 1.1920928955078125E-7);
        _delta_x[i] += _t121;
        _final_error += _t121;
    }
    clad::array<double> _delta_xexact(_d_xexact.size());
    i = 0;
    for (; i < _d_xexact.size(); i++) {
        double _t122 = std::abs(_d_xexact[i] * xexact[i] * 1.1920928955078125E-7);
        _delta_xexact[i] += _t122;
        _final_error += _t122;
    }
    i = 0;
    for (; i < _d_r.size(); i++) {
        double _t123 = std::abs(_d_r[i] * _EERepl_r0[i] * 1.1920928955078125E-7);
        _delta_r[i] += _t123;
        _final_error += _t123;
    }
    i = 0;
    for (; i < _d_p.size(); i++) {
        double _t124 = std::abs(_d_p[i] * _EERepl_p0[i] * 1.1920928955078125E-7);
        _delta_p[i] += _t124;
        _final_error += _t124;
    }
    clad::array<double> _delta_Ap(_d_Ap.size());
    i = 0;
    for (; i < _d_Ap.size(); i++) {
        double _t125 = std::abs(_d_Ap[i] * Ap[i] * 1.1920928955078125E-7);
        _delta_Ap[i] += _t125;
        _final_error += _t125;
    }
    _final_error += _delta_diff + _delta_oldrtrans + _delta_beta + _delta_sum + _delta_rtrans + _delta_tolerance + _delta_normr + _delta_alpha + _delta_residual;
}
