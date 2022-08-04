namespace clad {

HPC_Sparse_Matrix A;

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
    double _d_sum = 0;
    double _delta_sum = 0;
    double _EERepl_sum0;
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
    clad::tape<int> _t10 = {};
    clad::tape<double> _EERepl_sum1 = {};
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
    double _EERepl_rtrans1;
    unsigned long _t36;
    
    clad::tape<double> _t37 = {};
    clad::tape<int> _t38 = {};
    clad::tape<double> _t40 = {};
    clad::tape<int> _t41 = {};
    clad::tape<double> _EERepl_rtrans2 = {};
    double _t43;
    double _EERepl_normr2;
    unsigned long _t44;
    int _d_k = 0;
    clad::tape<bool> _t46 = {};
    clad::tape<double> _EERepl_beta2 = {};
    clad::tape<unsigned long> _t47 = {};
    
    clad::tape<int> _t48 = {};
    clad::tape<int> _t50 = {};
    clad::tape<double> _t52 = {};
    clad::tape<double> _t53 = {};
    clad::tape<int> _t54 = {};
    clad::tape<double> _EERepl_p2 = {};
    clad::tape<double> _EERepl_rtrans3 = {};
    clad::tape<unsigned long> _t56 = {};
    
    clad::tape<double> _t57 = {};
    clad::tape<int> _t58 = {};
    clad::tape<double> _t60 = {};
    clad::tape<int> _t61 = {};
    clad::tape<double> _EERepl_rtrans4 = {};
    clad::tape<double> _t63 = {};
    clad::tape<double> _t64 = {};
    clad::tape<double> _EERepl_beta3 = {};
    clad::tape<unsigned long> _t65 = {};
    
    clad::tape<int> _t66 = {};
    clad::tape<int> _t68 = {};
    clad::tape<double> _t70 = {};
    clad::tape<double> _t71 = {};
    clad::tape<int> _t72 = {};
    clad::tape<double> _EERepl_p3 = {};
    clad::tape<double> _t74 = {};
    clad::tape<double> _EERepl_normr3 = {};
    clad::tape<unsigned long> _t75 = {};
    clad::tape<int> _t76 = {};
    clad::tape<double> _EERepl_sum3 = {};
    clad::tape<unsigned long> _t78 = {};
    
    clad::tape<double> _t79 = {};
    clad::tape<int> _t80 = {};
    clad::tape<int> _t82 = {};
    clad::tape<double> _t84 = {};
    clad::tape<int> _t85 = {};
    clad::tape<int> _t87 = {};
    clad::tape<int> _t89 = {};
    clad::tape<double> _EERepl_sum4 = {};
    clad::tape<int> _t91 = {};
    clad::tape<double> _EERepl_beta4 = {};
    clad::tape<unsigned long> _t93 = {};
    
    clad::tape<double> _t94 = {};
    clad::tape<int> _t95 = {};
    clad::tape<double> _t97 = {};
    clad::tape<int> _t98 = {};
    clad::tape<double> _EERepl_beta5 = {};
    clad::tape<double> _t100 = {};
    clad::tape<double> _t101 = {};
    clad::tape<double> _EERepl_beta6 = {};
    clad::tape<unsigned long> _t102 = {};
    
    clad::tape<int> _t103 = {};
    clad::array<double> _delta_x(_d_x.size());
    clad::array<double> _EERepl_x0(_d_x.size());
    for (int i0 = 0; i0 < _d_x.size(); i0++) {
        _EERepl_x0[i0] = x[i0];
    }
    clad::tape<int> _t105 = {};
    clad::tape<double> _t107 = {};
    clad::tape<double> _t108 = {};
    clad::tape<int> _t109 = {};
    clad::tape<double> _EERepl_x1 = {};
    clad::tape<double> _EERepl_beta7 = {};
    clad::tape<unsigned long> _t111 = {};
    
    clad::tape<int> _t112 = {};
    clad::tape<int> _t114 = {};
    clad::tape<double> _t116 = {};
    clad::tape<double> _t117 = {};
    clad::tape<int> _t118 = {};
    clad::tape<double> _EERepl_r2 = {};
    double _d_residual = 0;
    double _delta_residual = 0;
    double _EERepl_residual0;
    double _d_diff = 0;
    double _delta_diff = 0;
    double _EERepl_diff0;
    unsigned long _t120;
    
    clad::tape<int> _t121 = {};
    clad::tape<int> _t123 = {};
    clad::tape<double> _t125 = {};
    clad::tape<double> _EERepl_diff1 = {};
    clad::tape<bool> _t128 = {};
    int niters = 0;
    double normr = 0.;
    _EERepl_normr0 = normr;
    int max_iter = 100;
    double tolerance = 0.;
    _EERepl_tolerance0 = tolerance;
    int cur_nnz;
    double sum;
    _EERepl_sum0 = sum;
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
        sum = 0;
        clad::push(_EERepl_sum1, sum);
        clad::push(_t12, 0UL);
        for (int j = 0; j < cur_nnz; j++) {
            clad::back(_t12)++;
            sum += clad::push(_t18, A.ptr_to_vals_in_row[clad::push(_t14, i)][clad::push(_t16, j)]) * clad::push(_t13, p[clad::push(_t23, A.ptr_to_inds_in_row[clad::push(_t19, i)][clad::push(_t21, j)])]);
            clad::push(_EERepl_sum2, sum);
        }
        Ap[clad::push(_t25, i)] = sum;
    }
    beta = -1.;
    _EERepl_beta1 = beta;
    _t27 = 0;
    for (int i = 0; i < nrow; i++) {
        _t27++;
        r[clad::push(_t28, i)] = b[clad::push(_t30, i)] + clad::push(_t33, beta) * clad::push(_t32, Ap[clad::push(_t34, i)]);
        clad::push(_EERepl_r1, r[clad::push(_t28, i)]);
    }
    rtrans = 0.;
    _EERepl_rtrans1 = rtrans;
    _t36 = 0;
    for (int i = 0; i < nrow; i++) {
        _t36++;
        rtrans += clad::push(_t40, r[clad::push(_t38, i)]) * clad::push(_t37, r[clad::push(_t41, i)]);
        clad::push(_EERepl_rtrans2, rtrans);
    }
    _t43 = rtrans;
    normr = sqrt(_t43);
    _EERepl_normr2 = normr;
    _t44 = 0;
    for (int k = 1; k < max_iter && normr > tolerance; k++) {
        _t44++;
        bool _t45 = k == 1;
        {
            if (_t45) {
                beta = 0.;
                clad::push(_EERepl_beta2, beta);
                clad::push(_t47, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t47)++;
                    p[clad::push(_t48, i)] = r[clad::push(_t50, i)] + clad::push(_t53, beta) * clad::push(_t52, r[clad::push(_t54, i)]);
                    clad::push(_EERepl_p2, p[clad::push(_t48, i)]);
                }
            } else {
                oldrtrans = rtrans;
                rtrans = 0.;
                clad::push(_EERepl_rtrans3, rtrans);
                clad::push(_t56, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t56)++;
                    rtrans += clad::push(_t60, r[clad::push(_t58, i)]) * clad::push(_t57, r[clad::push(_t61, i)]);
                    clad::push(_EERepl_rtrans4, rtrans);
                }
                beta = clad::push(_t64, rtrans) / clad::push(_t63, oldrtrans);
                clad::push(_EERepl_beta3, beta);
                clad::push(_t65, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t65)++;
                    p[clad::push(_t66, i)] = r[clad::push(_t68, i)] + clad::push(_t71, beta) * clad::push(_t70, p[clad::push(_t72, i)]);
                    clad::push(_EERepl_p3, p[clad::push(_t66, i)]);
                }
            }
            clad::push(_t46, _t45);
        }
        normr = sqrt(clad::push(_t74, rtrans));
        clad::push(_EERepl_normr3, normr);
        clad::push(_t75, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t75)++;
            cur_nnz = A.nnz_in_row[clad::push(_t76, i)];
            sum = 0;
            clad::push(_EERepl_sum3, sum);
            clad::push(_t78, 0UL);
            for (int j = 0; j < cur_nnz; j++) {
                clad::back(_t78)++;
                sum += clad::push(_t84, A.ptr_to_vals_in_row[clad::push(_t80, i)][clad::push(_t82, j)]) * clad::push(_t79, p[clad::push(_t89, A.ptr_to_inds_in_row[clad::push(_t85, i)][clad::push(_t87, j)])]);
                clad::push(_EERepl_sum4, sum);
            }
            Ap[clad::push(_t91, i)] = sum;
        }
        beta = 0.;
        clad::push(_EERepl_beta4, beta);
        clad::push(_t93, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t93)++;
            beta += clad::push(_t97, p[clad::push(_t95, i)]) * clad::push(_t94, Ap[clad::push(_t98, i)]);
            clad::push(_EERepl_beta5, beta);
        }
        beta = clad::push(_t101, rtrans) / clad::push(_t100, beta);
        clad::push(_EERepl_beta6, beta);
        clad::push(_t102, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t102)++;
            x[clad::push(_t103, i)] = x[clad::push(_t105, i)] + clad::push(_t108, beta) * clad::push(_t107, p[clad::push(_t109, i)]);
            clad::push(_EERepl_x1, x[clad::push(_t103, i)]);
        }
        beta = -beta;
        clad::push(_EERepl_beta7, beta);
        clad::push(_t111, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t111)++;
            r[clad::push(_t112, i)] = r[clad::push(_t114, i)] + clad::push(_t117, beta) * clad::push(_t116, Ap[clad::push(_t118, i)]);
            clad::push(_EERepl_r2, r[clad::push(_t112, i)]);
        }
        niters = k;
    }
    double residual = 0.;
    _EERepl_residual0 = residual;
    double diff = 0;
    _EERepl_diff0 = diff;
    _t120 = 0;
    for (int i = 0; i < nrow; i++) {
        _t120++;
        diff = fabs(clad::push(_t125, x[clad::push(_t121, i)] - xexact[clad::push(_t123, i)]));
        clad::push(_EERepl_diff1, diff);
        bool _t127 = diff > residual;
        {
            if (_t127)
                residual = diff;
            clad::push(_t128, _t127);
        }
    }
    double HPCCG_residual_return = residual;
    goto _label0;
  _label0:
    _d_residual += 1;
    for (; _t120; _t120--) {
        if (clad::pop(_t128)) {
            double _r_d31 = _d_residual;
            _d_diff += _r_d31;
            _d_residual -= _r_d31;
        }
        {
            double _r_d30 = _d_diff;
            double _r56 = _r_d30 * 1;
            int _t122 = clad::pop(_t121);
            _d_x[_t122] += _r56;
            int _t124 = clad::pop(_t123);
            _d_xexact[_t124] += -_r56;
            double _r57 = clad::pop(_EERepl_diff1);
            _delta_diff += clad::getErrorVal(_r_d30, _r57, "diff");
            _d_diff -= _r_d30;
        }
    }
    _delta_diff += clad::getErrorVal(_d_diff, _EERepl_diff0, "diff");
    _delta_residual += clad::getErrorVal(_d_residual, _EERepl_residual0, "residual");
    for (; _t44; _t44--) {
        {
            int _r_d29 = _d_niters;
            _d_k += _r_d29;
            _d_niters -= _r_d29;
        }
        {
            for (; clad::back(_t111); clad::back(_t111)--) {
                int _t113 = clad::pop(_t112);
                double _r_d28 = _d_r[_t113];
                int _t115 = clad::pop(_t114);
                _d_r[_t115] += _r_d28;
                double _r52 = _r_d28 * clad::pop(_t116);
                _d_beta += _r52;
                double _r53 = clad::pop(_t117) * _r_d28;
                int _t119 = clad::pop(_t118);
                _d_Ap[_t119] += _r53;
                double _r54 = clad::pop(_EERepl_r2);
                int _r55 = clad::pop(_t112);
                _delta_r[_r55] += clad::getErrorVal(_r_d28, _r54, "r");
                _final_error += _delta_r[_r55];
                _d_r[_t113] -= _r_d28;
                _d_r[_t113];
            }
            clad::pop(_t111);
        }
        {
            double _r_d27 = _d_beta;
            _d_beta += -_r_d27;
            double _r51 = clad::pop(_EERepl_beta7);
            _delta_beta += clad::getErrorVal(_r_d27, _r51, "beta");
            _d_beta -= _r_d27;
        }
        {
            for (; clad::back(_t102); clad::back(_t102)--) {
                int _t104 = clad::pop(_t103);
                double _r_d26 = _d_x[_t104];
                int _t106 = clad::pop(_t105);
                _d_x[_t106] += _r_d26;
                double _r47 = _r_d26 * clad::pop(_t107);
                _d_beta += _r47;
                double _r48 = clad::pop(_t108) * _r_d26;
                int _t110 = clad::pop(_t109);
                _d_p[_t110] += _r48;
                double _r49 = clad::pop(_EERepl_x1);
                int _r50 = clad::pop(_t103);
                _delta_x[_r50] += clad::getErrorVal(_r_d26, _r49, "x");
                _final_error += _delta_x[_r50];
                _d_x[_t104] -= _r_d26;
                _d_x[_t104];
            }
            clad::pop(_t102);
        }
        {
            double _r_d25 = _d_beta;
            double _r43 = clad::pop(_t100);
            double _r44 = _r_d25 / _r43;
            _d_rtrans += _r44;
            double _r45 = _r_d25 * -clad::pop(_t101) / (_r43 * _r43);
            _d_beta += _r45;
            double _r46 = clad::pop(_EERepl_beta6);
            _delta_beta += clad::getErrorVal(_r_d25, _r46, "beta");
            _d_beta -= _r_d25;
        }
        {
            for (; clad::back(_t93); clad::back(_t93)--) {
                double _r_d24 = _d_beta;
                _d_beta += _r_d24;
                double _r40 = _r_d24 * clad::pop(_t94);
                int _t96 = clad::pop(_t95);
                _d_p[_t96] += _r40;
                double _r41 = clad::pop(_t97) * _r_d24;
                int _t99 = clad::pop(_t98);
                _d_Ap[_t99] += _r41;
                double _r42 = clad::pop(_EERepl_beta5);
                _delta_beta += clad::getErrorVal(_r_d24, _r42, "beta");
                _d_beta -= _r_d24;
            }
            clad::pop(_t93);
        }
        {
            double _r_d23 = _d_beta;
            double _r39 = clad::pop(_EERepl_beta4);
            _delta_beta += clad::getErrorVal(_r_d23, _r39, "beta");
            _d_beta -= _r_d23;
        }
        {
            for (; clad::back(_t75); clad::back(_t75)--) {
                {
                    int _t92 = clad::pop(_t91);
                    int _t88 = clad::pop(_t87);
                    int _t86 = clad::pop(_t85);
                    double _r_d22 = _d_Ap[_t92];
                    _d_sum += _r_d22;
                    _d_Ap[_t92] -= _r_d22;
                    _d_Ap[_t92];
                }
                {
                    for (; clad::back(_t78); clad::back(_t78)--) {
                        double _r_d21 = _d_sum;
                        _d_sum += _r_d21;
                        double _r36 = _r_d21 * clad::pop(_t79);
                        int _t81 = clad::pop(_t80);
                        int _t83 = clad::pop(_t82);
                        double _r37 = clad::pop(_t84) * _r_d21;
                        int _t90 = clad::pop(_t89);
                        _d_p[_t90] += _r37;
                        double _r38 = clad::pop(_EERepl_sum4);
                        _delta_sum += clad::getErrorVal(_r_d21, _r38, "sum");
                        _d_sum -= _r_d21;
                    }
                    clad::pop(_t78);
                }
                {
                    double _r_d20 = _d_sum;
                    double _r35 = clad::pop(_EERepl_sum3);
                    _delta_sum += clad::getErrorVal(_r_d20, _r35, "sum");
                    _d_sum -= _r_d20;
                }
                {
                    int _r_d19 = _d_cur_nnz;
                    int _t77 = clad::pop(_t76);
                    _d_cur_nnz -= _r_d19;
                }
            }
            clad::pop(_t75);
        }
        {
            double _r_d18 = _d_normr;
            double _r33 = _r_d18 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t74), 1.).pushforward;
            _d_rtrans += _r33;
            double _r34 = clad::pop(_EERepl_normr3);
            _delta_normr += clad::getErrorVal(_r_d18, _r34, "normr");
            _d_normr -= _r_d18;
        }
        if (clad::pop(_t46)) {
            {
                for (; clad::back(_t47); clad::back(_t47)--) {
                    int _t49 = clad::pop(_t48);
                    double _r_d12 = _d_p[_t49];
                    int _t51 = clad::pop(_t50);
                    _d_r[_t51] += _r_d12;
                    double _r17 = _r_d12 * clad::pop(_t52);
                    _d_beta += _r17;
                    double _r18 = clad::pop(_t53) * _r_d12;
                    int _t55 = clad::pop(_t54);
                    _d_r[_t55] += _r18;
                    double _r19 = clad::pop(_EERepl_p2);
                    int _r20 = clad::pop(_t48);
                    _delta_p[_r20] += clad::getErrorVal(_r_d12, _r19, "p");
                    _final_error += _delta_p[_r20];
                    _d_p[_t49] -= _r_d12;
                    _d_p[_t49];
                }
                clad::pop(_t47);
            }
            {
                double _r_d11 = _d_beta;
                double _r16 = clad::pop(_EERepl_beta2);
                _delta_beta += clad::getErrorVal(_r_d11, _r16, "beta");
                _d_beta -= _r_d11;
            }
        } else {
            {
                for (; clad::back(_t65); clad::back(_t65)--) {
                    int _t67 = clad::pop(_t66);
                    double _r_d17 = _d_p[_t67];
                    int _t69 = clad::pop(_t68);
                    _d_r[_t69] += _r_d17;
                    double _r29 = _r_d17 * clad::pop(_t70);
                    _d_beta += _r29;
                    double _r30 = clad::pop(_t71) * _r_d17;
                    int _t73 = clad::pop(_t72);
                    _d_p[_t73] += _r30;
                    double _r31 = clad::pop(_EERepl_p3);
                    int _r32 = clad::pop(_t66);
                    _delta_p[_r32] += clad::getErrorVal(_r_d17, _r31, "p");
                    _final_error += _delta_p[_r32];
                    _d_p[_t67] -= _r_d17;
                    _d_p[_t67];
                }
                clad::pop(_t65);
            }
            {
                double _r_d16 = _d_beta;
                double _r25 = clad::pop(_t63);
                double _r26 = _r_d16 / _r25;
                _d_rtrans += _r26;
                double _r27 = _r_d16 * -clad::pop(_t64) / (_r25 * _r25);
                _d_oldrtrans += _r27;
                double _r28 = clad::pop(_EERepl_beta3);
                _delta_beta += clad::getErrorVal(_r_d16, _r28, "beta");
                _d_beta -= _r_d16;
            }
            {
                for (; clad::back(_t56); clad::back(_t56)--) {
                    double _r_d15 = _d_rtrans;
                    _d_rtrans += _r_d15;
                    double _r22 = _r_d15 * clad::pop(_t57);
                    int _t59 = clad::pop(_t58);
                    _d_r[_t59] += _r22;
                    double _r23 = clad::pop(_t60) * _r_d15;
                    int _t62 = clad::pop(_t61);
                    _d_r[_t62] += _r23;
                    double _r24 = clad::pop(_EERepl_rtrans4);
                    _delta_rtrans += clad::getErrorVal(_r_d15, _r24, "rtrans");
                    _d_rtrans -= _r_d15;
                }
                clad::pop(_t56);
            }
            {
                double _r_d14 = _d_rtrans;
                double _r21 = clad::pop(_EERepl_rtrans3);
                _delta_rtrans += clad::getErrorVal(_r_d14, _r21, "rtrans");
                _d_rtrans -= _r_d14;
            }
            {
                double _r_d13 = _d_oldrtrans;
                _d_rtrans += _r_d13;
                _d_oldrtrans -= _r_d13;
            }
        }
    }
    {
        double _r_d10 = _d_normr;
        double _r15 = _r_d10 * clad::custom_derivatives::sqrt_pushforward(_t43, 1.).pushforward;
        _d_rtrans += _r15;
        _delta_normr += clad::getErrorVal(_r_d10, _EERepl_normr2, "normr");
        _d_normr -= _r_d10;
    }
    for (; _t36; _t36--) {
        double _r_d9 = _d_rtrans;
        _d_rtrans += _r_d9;
        double _r12 = _r_d9 * clad::pop(_t37);
        int _t39 = clad::pop(_t38);
        _d_r[_t39] += _r12;
        double _r13 = clad::pop(_t40) * _r_d9;
        int _t42 = clad::pop(_t41);
        _d_r[_t42] += _r13;
        double _r14 = clad::pop(_EERepl_rtrans2);
        _delta_rtrans += clad::getErrorVal(_r_d9, _r14, "rtrans");
        _d_rtrans -= _r_d9;
    }
    {
        double _r_d8 = _d_rtrans;
        _delta_rtrans += clad::getErrorVal(_r_d8, _EERepl_rtrans1, "rtrans");
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
        _delta_r[_r11] += clad::getErrorVal(_r_d7, _r10, "r");
        _final_error += _delta_r[_r11];
        _d_r[_t29] -= _r_d7;
        _d_r[_t29];
    }
    {
        double _r_d6 = _d_beta;
        _delta_beta += clad::getErrorVal(_r_d6, _EERepl_beta1, "beta");
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
                _delta_sum += clad::getErrorVal(_r_d4, _r7, "sum");
                _d_sum -= _r_d4;
            }
            clad::pop(_t12);
        }
        {
            double _r_d3 = _d_sum;
            double _r4 = clad::pop(_EERepl_sum1);
            _delta_sum += clad::getErrorVal(_r_d3, _r4, "sum");
            _d_sum -= _r_d3;
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
        double _t129 = clad::getErrorVal(_d_b[i], b[i], "b");
        _delta_b[i] += _t129;
        _final_error += _t129;
    }
    i = 0;
    double errSum = 0;
    for (; i < _d_x.size(); i++) {
        double _t130 = clad::getErrorVal(_d_x[i], _EERepl_x0[i], "x");
        _delta_x[i] += _t130;
        errSum += _delta_x[i];
        _final_error += _t130;
    }
    clad::array<double> _delta_xexact(_d_xexact.size());
    i = 0;
    for (; i < _d_xexact.size(); i++) {
        double _t131 = clad::getErrorVal(_d_xexact[i], xexact[i], "xexact");
        _delta_xexact[i] += _t131;
        _final_error += _t131;
    }
    i = 0;
    double errSumr = 0;
    for (; i < _d_r.size(); i++) {
        double _t132 = clad::getErrorVal(_d_r[i], _EERepl_r0[i], "r");
        _delta_r[i] += _t132;
        errSumr += _delta_r[i];
        _final_error += _t132;
    }
    i = 0;
    for (; i < _d_p.size(); i++) {
        double _t133 = clad::getErrorVal(_d_p[i], _EERepl_p0[i], "p");
        _delta_p[i] += _t133;
        _final_error += _t133;
    }
    clad::array<double> _delta_Ap(_d_Ap.size());
    i = 0;
    for (; i < _d_Ap.size(); i++) {
        double _t134 = clad::getErrorVal(_d_Ap[i], Ap[i], "Ap");
        _delta_Ap[i] += _t134;
        _final_error += _t134;
    }

    std::cout << "sum : " << _delta_sum << std::endl;
    std::cout << "x : " << errSum << std::endl;
    std::cout << "r : " << errSumr << std::endl;
    
    _final_error += _delta_normr + _delta_sum + _delta_rtrans + _delta_beta + _delta_diff + _delta_tolerance + _delta_residual + _delta_oldrtrans;
}

}

