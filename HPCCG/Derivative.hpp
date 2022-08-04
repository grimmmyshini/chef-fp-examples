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
    double _d_sum = 0, _d_temp = 0;
    double _delta_sum = 0;
    double _EERepl_sum0;
    double _delta_sum0 = 0;
    double _EERepl_sum1;
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
    clad::tape<int> _t3 = {};
    clad::array<double> _delta_p(_d_p.size());
    clad::array<double> _EERepl_p0(_d_p.size());
    for (int i0 = 0; i0 < _d_p.size(); i0++) {
        _EERepl_p0[i0] = p[i0];
    }
    clad::tape<double> _t5 = {};
    clad::tape<double> _t6 = {};
    clad::tape<double> _EERepl_p1 = {};
    unsigned long _t7;
    
    clad::tape<int> _t8 = {};
    clad::tape<double> _EERepl_sum2 = {};
    clad::tape<unsigned long> _t10 = {};

    clad::tape<double> _t11 = {};
    clad::tape<int> _t12 = {};
    clad::tape<int> _t14 = {};
    clad::tape<double> _t16 = {};
    clad::tape<int> _t17 = {};
    clad::tape<int> _t19 = {};
    clad::tape<int> _t21 = {};
    clad::tape<double> _EERepl_sum3 = {};
    clad::tape<int> _t23 = {};
    double _EERepl_beta1;
    unsigned long _t25;
    
    clad::tape<int> _t26 = {};
    clad::array<double> _delta_r(_d_r.size());
    clad::array<double> _EERepl_r0(_d_r.size());
    for (int i0 = 0; i0 < _d_r.size(); i0++) {
        _EERepl_r0[i0] = r[i0];
    }
    clad::tape<int> _t28 = {};
    clad::tape<double> _t30 = {};
    clad::tape<double> _t31 = {};
    clad::tape<int> _t32 = {};
    clad::tape<double> _EERepl_r1 = {};
    double _EERepl_rtrans1;
    unsigned long _t34;
    
    clad::tape<int> _t35 = {};
    clad::tape<double> _t37 = {};
    clad::tape<double> _t38 = {};
    clad::tape<double> _EERepl_rtrans2 = {};
    double _t39;
    double _EERepl_normr2;
    unsigned long _t40;
    int _d_k = 0;
    clad::tape<bool> _t42 = {};
    clad::tape<double> _EERepl_beta2 = {};
    clad::tape<unsigned long> _t43 = {};
    
    clad::tape<int> _t44 = {};
    clad::tape<int> _t46 = {};
    clad::tape<double> _t48 = {};
    clad::tape<double> _t49 = {};
    clad::tape<double> _EERepl_p2 = {};
    clad::tape<double> _EERepl_rtrans3 = {};
    clad::tape<unsigned long> _t50 = {};
    
    clad::tape<int> _t51 = {};
    clad::tape<double> _t53 = {};
    clad::tape<double> _t54 = {};
    clad::tape<double> _EERepl_rtrans4 = {};
    clad::tape<double> _t55 = {};
    clad::tape<double> _t56 = {};
    clad::tape<double> _EERepl_beta3 = {};
    clad::tape<unsigned long> _t57 = {};
    
    clad::tape<int> _t58 = {};
    clad::tape<int> _t60 = {};
    clad::tape<double> _t62 = {};
    clad::tape<double> _t63 = {};
    clad::tape<int> _t64 = {};
    clad::tape<double> _EERepl_p3 = {};
    clad::tape<double> _t66 = {};
    clad::tape<double> _EERepl_normr3 = {};
    clad::tape<unsigned long> _t67 = {};
    
    clad::tape<int> _t68 = {};
    clad::tape<double> _EERepl_sum4 = {};
    clad::tape<unsigned long> _t70 = {};

    clad::tape<double> _t71 = {};
    clad::tape<int> _t72 = {};
    clad::tape<int> _t74 = {};
    clad::tape<double> _t76 = {};
    clad::tape<int> _t77 = {};
    clad::tape<int> _t79 = {};
    clad::tape<int> _t81 = {};
    clad::tape<double> _EERepl_sum5 = {};
    clad::tape<int> _t83 = {};
    clad::tape<double> _EERepl_beta4 = {};
    clad::tape<unsigned long> _t85 = {};
    
    clad::tape<double> _t86 = {};
    clad::tape<int> _t87 = {};
    clad::tape<double> _t89 = {};
    clad::tape<int> _t90 = {};
    clad::tape<double> _EERepl_beta5 = {};
    clad::tape<double> _t92 = {};
    clad::tape<double> _t93 = {};
    clad::tape<double> _EERepl_beta6 = {};
    clad::tape<unsigned long> _t94 = {};
    
    clad::tape<int> _t95 = {};
    clad::array<double> _delta_x(_d_x.size());
    clad::array<double> _EERepl_x0(_d_x.size());
    for (int i0 = 0; i0 < _d_x.size(); i0++) {
        _EERepl_x0[i0] = x[i0];
    }
    clad::tape<int> _t97 = {};
    clad::tape<double> _t99 = {};
    clad::tape<double> _t100 = {};
    clad::tape<int> _t101 = {};
    clad::tape<double> _EERepl_x1 = {};
    clad::tape<double> _EERepl_beta7 = {};
    clad::tape<unsigned long> _t103 = {};
    
    clad::tape<int> _t104 = {};
    clad::tape<int> _t106 = {};
    clad::tape<double> _t108 = {};
    clad::tape<double> _t109 = {};
    clad::tape<int> _t110 = {};
    clad::tape<double> _EERepl_r2 = {};
    double _d_residual = 0;
    double _delta_residual = 0;
    double _EERepl_residual0;
    double _d_diff = 0;
    double _delta_diff = 0;
    double _EERepl_diff0;
    unsigned long _t112;
    
    clad::tape<int> _t113 = {};
    clad::tape<int> _t115 = {};
    clad::tape<double> _t117 = {};
    clad::tape<double> _EERepl_diff1 = {};
    clad::tape<bool> _t120 = {};
    int niters = 0;
    double normr = 0.;
    _EERepl_normr0 = normr;
    int max_iter = 100;
    double tolerance = 0.;
    _EERepl_tolerance0 = tolerance;
    int cur_nnz;
    double sum, temp;
    _EERepl_sum1 = sum;
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
    // int _EERepl_cur_nnz1 = cur_nnz;
    // int _EERepl_cur_nnz0 = cur_nnz;
    for (int i = 0; i < nrow; i++) {
        _t0++;
        temp = x[clad::push(_t1, i)];
        p[clad::push(_t3, i)] = temp + clad::push(_t6, beta) * clad::push(_t5, temp);
        clad::push(_EERepl_p1, p[clad::push(_t3, i)]);
    }
    _t7 = 0;
    for (int i = 0; i < nrow; i++) {
        _t7++;
        cur_nnz = A.nnz_in_row[clad::push(_t8, i)];
        sum = 0;
        clad::push(_EERepl_sum2, sum);
        clad::push(_t10, 0UL);
        for (int j = 0; j < cur_nnz; j++) {
            clad::back(_t10)++;
            sum += clad::push(_t16, A.ptr_to_vals_in_row[clad::push(_t12, i)][clad::push(_t14, j)]) * clad::push(_t11, p[clad::push(_t21, A.ptr_to_inds_in_row[clad::push(_t17, i)][clad::push(_t19, j)])]);
            clad::push(_EERepl_sum3, sum);
        }
        Ap[clad::push(_t23, i)] = sum;
    }
    beta = -1.;
    _EERepl_beta1 = beta;
    _t25 = 0;
    for (int i = 0; i < nrow; i++) {
        _t25++;
        r[clad::push(_t26, i)] = b[clad::push(_t28, i)] + clad::push(_t31, beta) * clad::push(_t30, Ap[clad::push(_t32, i)]);
        clad::push(_EERepl_r1, r[clad::push(_t26, i)]);
    }
    rtrans = 0.;
    _EERepl_rtrans1 = rtrans;
    _t34 = 0;
    for (int i = 0; i < nrow; i++) {
        _t34++;
        temp = r[clad::push(_t35, i)];
        rtrans += clad::push(_t38, temp) * clad::push(_t37, temp);
        clad::push(_EERepl_rtrans2, rtrans);
    }
    _t39 = rtrans;
    normr = sqrt(_t39);
    _EERepl_normr2 = normr;
    _t40 = 0;
    for (int k = 1; k < max_iter && normr > tolerance; k++) {
        _t40++;
        bool _t41 = k == 1;
        {
            if (_t41) {
                beta = 0.;
                clad::push(_EERepl_beta2, beta);
                clad::push(_t43, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t43)++;
                    temp = r[clad::push(_t44, i)];
                    p[clad::push(_t46, i)] = temp + clad::push(_t49, beta) * clad::push(_t48, temp);
                    clad::push(_EERepl_p2, p[clad::push(_t46, i)]);
                }
            } else {
                oldrtrans = rtrans;
                rtrans = 0.;
                clad::push(_EERepl_rtrans3, rtrans);
                clad::push(_t50, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t50)++;
                    temp = r[clad::push(_t51, i)];
                    rtrans += clad::push(_t54, temp) * clad::push(_t53, temp);
                    clad::push(_EERepl_rtrans4, rtrans);
                }
                beta = clad::push(_t56, rtrans) / clad::push(_t55, oldrtrans);
                clad::push(_EERepl_beta3, beta);
                clad::push(_t57, 0UL);
                for (int i = 0; i < nrow; i++) {
                    clad::back(_t57)++;
                    p[clad::push(_t58, i)] = r[clad::push(_t60, i)] + clad::push(_t63, beta) * clad::push(_t62, p[clad::push(_t64, i)]);
                    clad::push(_EERepl_p3, p[clad::push(_t58, i)]);
                }
            }
            clad::push(_t42, _t41);
        }
        normr = sqrt(clad::push(_t66, rtrans));
        clad::push(_EERepl_normr3, normr);
        clad::push(_t67, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t67)++;
            cur_nnz = A.nnz_in_row[clad::push(_t68, i)];
            sum = 0;
            clad::push(_EERepl_sum4, sum);
            clad::push(_t70, 0UL);
            for (int j = 0; j < cur_nnz; j++) {
                clad::back(_t70)++;
                sum += clad::push(_t76, A.ptr_to_vals_in_row[clad::push(_t72, i)][clad::push(_t74, j)]) * clad::push(_t71, p[clad::push(_t81, A.ptr_to_inds_in_row[clad::push(_t77, i)][clad::push(_t79, j)])]);
                clad::push(_EERepl_sum5, sum);
            }
            Ap[clad::push(_t83, i)] = sum;
        }
        beta = 0.;
        clad::push(_EERepl_beta4, beta);
        clad::push(_t85, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t85)++;
            beta += clad::push(_t89, p[clad::push(_t87, i)]) * clad::push(_t86, Ap[clad::push(_t90, i)]);
            clad::push(_EERepl_beta5, beta);
        }
        beta = clad::push(_t93, rtrans) / clad::push(_t92, beta);
        clad::push(_EERepl_beta6, beta);
        clad::push(_t94, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t94)++;
            x[clad::push(_t95, i)] = x[clad::push(_t97, i)] + clad::push(_t100, beta) * clad::push(_t99, p[clad::push(_t101, i)]);
            clad::push(_EERepl_x1, x[clad::push(_t95, i)]);
        }
        beta = -beta;
        clad::push(_EERepl_beta7, beta);
        clad::push(_t103, 0UL);
        for (int i = 0; i < nrow; i++) {
            clad::back(_t103)++;
            r[clad::push(_t104, i)] = r[clad::push(_t106, i)] + clad::push(_t109, beta) * clad::push(_t108, Ap[clad::push(_t110, i)]);
            clad::push(_EERepl_r2, r[clad::push(_t104, i)]);
        }
        niters = k;
    }
    double residual = 0.;
    _EERepl_residual0 = residual;
    double diff = 0;
    _EERepl_diff0 = diff;
    _t112 = 0;
    for (int i = 0; i < nrow; i++) {
        _t112++;
        diff = fabs(clad::push(_t117, x[clad::push(_t113, i)] - xexact[clad::push(_t115, i)]));
        clad::push(_EERepl_diff1, diff);
        bool _t119 = diff > residual;
        {
            if (_t119)
                residual = diff;
            clad::push(_t120, _t119);
        }
    }
    double HPCCG_residual_return = residual;
    goto _label0;
  _label0:
    _d_residual += 1;
    for (; _t112; _t112--) {
        if (clad::pop(_t120)) {
            double _r_d35 = _d_residual;
            _d_diff += _r_d35;
            _d_residual -= _r_d35;
        }
        {
            double _r_d34 = _d_diff;
            double _r56 = _r_d34 * 1;
            int _t114 = clad::pop(_t113);
            _d_x[_t114] += _r56;
            int _t116 = clad::pop(_t115);
            _d_xexact[_t116] += -_r56;
            double _r57 = clad::pop(_EERepl_diff1);
            _delta_diff += clad::getErrorVal(_r_d34, _r57, "diff");
            _d_diff -= _r_d34;
        }
    }
    _delta_diff += clad::getErrorVal(_d_diff, _EERepl_diff0, "diff");
    _delta_residual += clad::getErrorVal(_d_residual, _EERepl_residual0, "residual");
    for (; _t40; _t40--) {
        {
            int _r_d33 = _d_niters;
            _d_k += _r_d33;
            _d_niters -= _r_d33;
        }
        {
            for (; clad::back(_t103); clad::back(_t103)--) {
                int _t105 = clad::pop(_t104);
                double _r_d32 = _d_r[_t105];
                int _t107 = clad::pop(_t106);
                _d_r[_t107] += _r_d32;
                double _r52 = _r_d32 * clad::pop(_t108);
                _d_beta += _r52;
                double _r53 = clad::pop(_t109) * _r_d32;
                int _t111 = clad::pop(_t110);
                _d_Ap[_t111] += _r53;
                double _r54 = clad::pop(_EERepl_r2);
                int _r55 = clad::pop(_t104);
                _delta_r[_r55] += clad::getErrorVal(_r_d32, _r54, "r");
                _final_error += _delta_r[_r55];
                _d_r[_t105] -= _r_d32;
                _d_r[_t105];
            }
            clad::pop(_t103);
        }
        {
            double _r_d31 = _d_beta;
            _d_beta += -_r_d31;
            double _r51 = clad::pop(_EERepl_beta7);
            _delta_beta += clad::getErrorVal(_r_d31, _r51, "beta");
            _d_beta -= _r_d31;
        }
        {
            for (; clad::back(_t94); clad::back(_t94)--) {
                int _t96 = clad::pop(_t95);
                double _r_d30 = _d_x[_t96];
                int _t98 = clad::pop(_t97);
                _d_x[_t98] += _r_d30;
                double _r47 = _r_d30 * clad::pop(_t99);
                _d_beta += _r47;
                double _r48 = clad::pop(_t100) * _r_d30;
                int _t102 = clad::pop(_t101);
                _d_p[_t102] += _r48;
                double _r49 = clad::pop(_EERepl_x1);
                int _r50 = clad::pop(_t95);
                _delta_x[_r50] += clad::getErrorVal(_r_d30, _r49, "x");
                _final_error += _delta_x[_r50];
                _d_x[_t96] -= _r_d30;
                _d_x[_t96];
            }
            clad::pop(_t94);
        }
        {
            double _r_d29 = _d_beta;
            double _r43 = clad::pop(_t92);
            double _r44 = _r_d29 / _r43;
            _d_rtrans += _r44;
            double _r45 = _r_d29 * -clad::pop(_t93) / (_r43 * _r43);
            _d_beta += _r45;
            double _r46 = clad::pop(_EERepl_beta6);
            _delta_beta += clad::getErrorVal(_r_d29, _r46, "beta");
            _d_beta -= _r_d29;
        }
        {
            for (; clad::back(_t85); clad::back(_t85)--) {
                double _r_d28 = _d_beta;
                _d_beta += _r_d28;
                double _r40 = _r_d28 * clad::pop(_t86);
                int _t88 = clad::pop(_t87);
                _d_p[_t88] += _r40;
                double _r41 = clad::pop(_t89) * _r_d28;
                int _t91 = clad::pop(_t90);
                _d_Ap[_t91] += _r41;
                double _r42 = clad::pop(_EERepl_beta5);
                _delta_beta += clad::getErrorVal(_r_d28, _r42, "beta");
                _d_beta -= _r_d28;
            }
            clad::pop(_t85);
        }
        {
            double _r_d27 = _d_beta;
            double _r39 = clad::pop(_EERepl_beta4);
            _delta_beta += clad::getErrorVal(_r_d27, _r39, "beta");
            _d_beta -= _r_d27;
        }
        {
            for (; clad::back(_t67); clad::back(_t67)--) {
                {
                    int _t84 = clad::pop(_t83);
                    int _t80 = clad::pop(_t79);
                    int _t78 = clad::pop(_t77);
                    double _r_d26 = _d_Ap[_t84];
                    _d_sum += _r_d26;
                    _d_Ap[_t84] -= _r_d26;
                    _d_Ap[_t84];
                }
                {
                    for (; clad::back(_t70); clad::back(_t70)--) {
                        double _r_d25 = _d_sum;
                        _d_sum += _r_d25;
                        double _r36 = _r_d25 * clad::pop(_t71);
                        int _t73 = clad::pop(_t72);
                        int _t75 = clad::pop(_t74);
                        double _r37 = clad::pop(_t76) * _r_d25;
                        int _t82 = clad::pop(_t81);
                        _d_p[_t82] += _r37;
                        double _r38 = clad::pop(_EERepl_sum5);
                        _delta_sum += clad::getErrorVal(_r_d25, _r38, "sum");
                        _d_sum -= _r_d25;
                    }
                    clad::pop(_t70);
                }
                {
                    double _r_d24 = _d_sum;
                    double _r35 = clad::pop(_EERepl_sum4);
                    _delta_sum += clad::getErrorVal(_r_d24, _r35, "sum");
                    _d_sum -= _r_d24;
                }
                {
                    int _r_d23 = _d_cur_nnz;
                    int _t69 = clad::pop(_t68);
                    _d_cur_nnz -= _r_d23;
                }
            }
            clad::pop(_t67);
        }
        {
            double _r_d22 = _d_normr;
            double _r33 = _r_d22 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t66), 1.).pushforward;
            _d_rtrans += _r33;
            double _r34 = clad::pop(_EERepl_normr3);
            _delta_normr += clad::getErrorVal(_r_d22, _r34, "normr");
            _d_normr -= _r_d22;
        }
        if (clad::pop(_t42)) {
            {
                for (; clad::back(_t43); clad::back(_t43)--) {
                    {
                        int _t47 = clad::pop(_t46);
                        double _r_d15 = _d_p[_t47];
                        _d_temp += _r_d15;
                        double _r17 = _r_d15 * clad::pop(_t48);
                        _d_beta += _r17;
                        double _r18 = clad::pop(_t49) * _r_d15;
                        _d_temp += _r18;
                        double _r19 = clad::pop(_EERepl_p2);
                        int _r20 = clad::pop(_t46);
                        _delta_p[_r20] += clad::getErrorVal(_r_d15, _r19, "p");
                        _final_error += _delta_p[_r20];
                        _d_p[_t47] -= _r_d15;
                        _d_p[_t47];
                    }
                    {
                        double _r_d14 = _d_temp;
                        int _t45 = clad::pop(_t44);
                        _d_r[_t45] += _r_d14;
                        _d_temp -= _r_d14;
                    }
                }
                clad::pop(_t43);
            }
            {
                double _r_d13 = _d_beta;
                double _r16 = clad::pop(_EERepl_beta2);
                _delta_beta += clad::getErrorVal(_r_d13, _r16, "beta");
                _d_beta -= _r_d13;
            }
        } else {
            {
                for (; clad::back(_t57); clad::back(_t57)--) {
                    int _t59 = clad::pop(_t58);
                    double _r_d21 = _d_p[_t59];
                    int _t61 = clad::pop(_t60);
                    _d_r[_t61] += _r_d21;
                    double _r29 = _r_d21 * clad::pop(_t62);
                    _d_beta += _r29;
                    double _r30 = clad::pop(_t63) * _r_d21;
                    int _t65 = clad::pop(_t64);
                    _d_p[_t65] += _r30;
                    double _r31 = clad::pop(_EERepl_p3);
                    int _r32 = clad::pop(_t58);
                    _delta_p[_r32] += clad::getErrorVal(_r_d21, _r31, "p");
                    _final_error += _delta_p[_r32];
                    _d_p[_t59] -= _r_d21;
                    _d_p[_t59];
                }
                clad::pop(_t57);
            }
            {
                double _r_d20 = _d_beta;
                double _r25 = clad::pop(_t55);
                double _r26 = _r_d20 / _r25;
                _d_rtrans += _r26;
                double _r27 = _r_d20 * -clad::pop(_t56) / (_r25 * _r25);
                _d_oldrtrans += _r27;
                double _r28 = clad::pop(_EERepl_beta3);
                _delta_beta += clad::getErrorVal(_r_d20, _r28, "beta");
                _d_beta -= _r_d20;
            }
            {
                for (; clad::back(_t50); clad::back(_t50)--) {
                    {
                        double _r_d19 = _d_rtrans;
                        _d_rtrans += _r_d19;
                        double _r22 = _r_d19 * clad::pop(_t53);
                        _d_temp += _r22;
                        double _r23 = clad::pop(_t54) * _r_d19;
                        _d_temp += _r23;
                        double _r24 = clad::pop(_EERepl_rtrans4);
                        _delta_rtrans += clad::getErrorVal(_r_d19, _r24, "rtrans");
                        _d_rtrans -= _r_d19;
                    }
                    {
                        double _r_d18 = _d_temp;
                        int _t52 = clad::pop(_t51);
                        _d_r[_t52] += _r_d18;
                        _d_temp -= _r_d18;
                    }
                }
                clad::pop(_t50);
            }
            {
                double _r_d17 = _d_rtrans;
                double _r21 = clad::pop(_EERepl_rtrans3);
                _delta_rtrans += clad::getErrorVal(_r_d17, _r21, "rtrans");
                _d_rtrans -= _r_d17;
            }
            {
                double _r_d16 = _d_oldrtrans;
                _d_rtrans += _r_d16;
                _d_oldrtrans -= _r_d16;
            }
        }
    }
    {
        double _r_d12 = _d_normr;
        double _r15 = _r_d12 * clad::custom_derivatives::sqrt_pushforward(_t39, 1.).pushforward;
        _d_rtrans += _r15;
        _delta_normr += clad::getErrorVal(_r_d12, _EERepl_normr2, "normr");
        _d_normr -= _r_d12;
    }
    for (; _t34; _t34--) {
        {
            double _r_d11 = _d_rtrans;
            _d_rtrans += _r_d11;
            double _r12 = _r_d11 * clad::pop(_t37);
            _d_temp += _r12;
            double _r13 = clad::pop(_t38) * _r_d11;
            _d_temp += _r13;
            double _r14 = clad::pop(_EERepl_rtrans2);
            _delta_rtrans += clad::getErrorVal(_r_d11, _r14, "rtrans");
            _d_rtrans -= _r_d11;
        }
        {
            double _r_d10 = _d_temp;
            int _t36 = clad::pop(_t35);
            _d_r[_t36] += _r_d10;
            _d_temp -= _r_d10;
        }
    }
    {
        double _r_d9 = _d_rtrans;
        _delta_rtrans += clad::getErrorVal(_r_d9, _EERepl_rtrans1, "rtrans");
        _d_rtrans -= _r_d9;
    }
    for (; _t25; _t25--) {
        int _t27 = clad::pop(_t26);
        double _r_d8 = _d_r[_t27];
        int _t29 = clad::pop(_t28);
        _d_b[_t29] += _r_d8;
        double _r8 = _r_d8 * clad::pop(_t30);
        _d_beta += _r8;
        double _r9 = clad::pop(_t31) * _r_d8;
        int _t33 = clad::pop(_t32);
        _d_Ap[_t33] += _r9;
        double _r10 = clad::pop(_EERepl_r1);
        int _r11 = clad::pop(_t26);
        _delta_r[_r11] += clad::getErrorVal(_r_d8, _r10, "r");
        _final_error += _delta_r[_r11];
        _d_r[_t27] -= _r_d8;
        _d_r[_t27];
    }
    {
        double _r_d7 = _d_beta;
        _delta_beta += clad::getErrorVal(_r_d7, _EERepl_beta1, "beta");
        _d_beta -= _r_d7;
    }
    for (; _t7; _t7--) {
        {
            int _t24 = clad::pop(_t23);
            int _t20 = clad::pop(_t19);
            int _t18 = clad::pop(_t17);
            double _r_d6 = _d_Ap[_t24];
            _d_sum += _r_d6;
            _d_Ap[_t24] -= _r_d6;
            _d_Ap[_t24];
        }
        {
            for (; clad::back(_t10); clad::back(_t10)--) {
                double _r_d5 = _d_sum;
                _d_sum += _r_d5;
                double _r5 = _r_d5 * clad::pop(_t11);
                int _t13 = clad::pop(_t12);
                int _t15 = clad::pop(_t14);
                double _r6 = clad::pop(_t16) * _r_d5;
                int _t22 = clad::pop(_t21);
                _d_p[_t22] += _r6;
                double _r7 = clad::pop(_EERepl_sum3);
                _delta_sum += clad::getErrorVal(_r_d5, _r7, "sum");
                _d_sum -= _r_d5;
            }
            clad::pop(_t10);
        }
        {
            double _r_d4 = _d_sum;
            double _r4 = clad::pop(_EERepl_sum2);
            _delta_sum += clad::getErrorVal(_r_d4, _r4, "sum");
            _d_sum -= _r_d4;
        }
        {
            int _r_d3 = _d_cur_nnz;
            int _t9 = clad::pop(_t8);
            _d_cur_nnz -= _r_d3;
        }
    }
    for (; _t0; _t0--) {
            int _t4 = clad::pop(_t3);
            double _r_d2 = _d_p[_t4];
            _d_temp += _r_d2;
            double _r0 = _r_d2 * clad::pop(_t5);
            _d_beta += _r0;
            double _r1 = clad::pop(_t6) * _r_d2;
            _d_temp += _r1;
            double _r2 = clad::pop(_EERepl_p1);
            int _r3 = clad::pop(_t3);
            _delta_p[_r3] += clad::getErrorVal(_r_d2, _r2, "p");
            _final_error += _delta_p[_r3];
            _d_p[_t4] -= _r_d2;
            _d_p[_t4];
            double _r_d1 = _d_temp;
            int _t2 = clad::pop(_t1);
            _d_x[_t2] += _r_d1;
            _d_temp -= _r_d1;
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
        double _t121 = clad::getErrorVal(_d_b[i], b[i], "b");
        _delta_b[i] += _t121;
        _final_error += _t121;
    }
    i = 0;
    double errSum = 0;
    for (; i < _d_x.size(); i++) {
        double _t122 = clad::getErrorVal(_d_x[i], _EERepl_x0[i], "x");
        _delta_x[i] += _t122;
        errSum += _delta_x[i];
        _final_error += _t122;
    }
    clad::array<double> _delta_xexact(_d_xexact.size());
    i = 0;
    for (; i < _d_xexact.size(); i++) {
        double _t123 = clad::getErrorVal(_d_xexact[i], xexact[i], "xexact");
        _delta_xexact[i] += _t123;
        _final_error += _t123;
    }
    i = 0;
    double errSumr = 0;
    for (; i < _d_r.size(); i++) {
        double _t124 = clad::getErrorVal(_d_r[i], _EERepl_r0[i], "r");
        _delta_r[i] += _t124;
        errSumr += _delta_r[i];
        _final_error += _t124;
    }
    i = 0;
    for (; i < _d_p.size(); i++) {
        double _t125 = clad::getErrorVal(_d_p[i], _EERepl_p0[i], "p");
        _delta_p[i] += _t125;
        _final_error += _t125;
    }
    clad::array<double> _delta_Ap(_d_Ap.size());
    i = 0;
    for (; i < _d_Ap.size(); i++) {
        double _t126 = clad::getErrorVal(_d_Ap[i], Ap[i], "Ap");
        _delta_Ap[i] += _t126;
        _final_error += _t126;
    }
    std::cout << "sum : " << _delta_sum << std::endl;
    std::cout << "x : " << errSum << std::endl;
    std::cout << "r : " << errSumr << std::endl;
    

    _final_error += _delta_normr + _delta_tolerance + _delta_rtrans + _delta_oldrtrans + _delta_sum + _delta_beta + _delta_residual + _delta_diff;
}

}