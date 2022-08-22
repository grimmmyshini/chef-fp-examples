namespace clad
{

    HPC_Sparse_Matrix A;

    void HPCCG_grad(double *b, double *x, int max_iter, double tolerance,
                    int &niters, double &normr, double *r, double *p,
                    double *Ap, double *xexact, clad::array_ref<double> _d_b,
                    clad::array_ref<double> _d_x, clad::array_ref<int> _d_max_iter,
                    clad::array_ref<double> _d_tolerance, clad::array_ref<int> _d_niters,
                    clad::array_ref<double> _d_normr, clad::array_ref<double> _d_r,
                    clad::array_ref<double> _d_p, clad::array_ref<double> _d_Ap,
                    clad::array_ref<double> _d_xexact, double &_final_error)
    {
        int _d_nrow = 0;
        int _d_ncol = 0;
        double _d_rtrans = 0;
        double _delta_rtrans = 0;
        double _EERepl_rtrans0;
        double _d_oldrtrans = 0;
        double _delta_oldrtrans = 0;
        double _EERepl_oldrtrans0;
        double _d_alp = 0, _d_bta = 0;
        double _delta_alp = 0;
        double _EERepl_alp0;
        double _delta_alp0 = 0;
        double _EERepl_alp1;
        int _d_rank = 0;
        int _d_print_freq = 0;
        bool _cond0;
        bool _cond1;
        double _EERepl_alp2;
        double _delta_bta = 0;
        double _EERepl_bta1;
        bool _cond2;
        unsigned long _t0;
        clad::tape<int> _t1 = {};
        clad::array<double> _delta_p(_d_p.size());
        clad::array<double> _EERepl_p0(_d_p.size());
        for (int i0 = 0; i0 < _d_p.size(); i0++)
        {
            _EERepl_p0[i0] = p[i0];
        }
        clad::tape<int> _t3 = {};
        clad::tape<double> _t5 = {};
        clad::tape<double> _t6 = {};
        clad::tape<int> _t7 = {};
        clad::tape<double> _EERepl_p1 = {};
        bool _cond3;
        unsigned long _t9;
        clad::tape<int> _t10 = {};
        clad::tape<double> _t12 = {};
        clad::tape<double> _t13 = {};
        clad::tape<int> _t14 = {};
        clad::tape<int> _t16 = {};
        clad::tape<double> _EERepl_p2 = {};
        unsigned long _t18;
        clad::tape<int> _t19 = {};
        clad::tape<double> _t21 = {};
        clad::tape<double> _t22 = {};
        clad::tape<int> _t23 = {};
        clad::tape<double> _t25 = {};
        clad::tape<double> _t26 = {};
        clad::tape<int> _t27 = {};
        clad::tape<double> _EERepl_p3 = {};
        unsigned long _t29;
        clad::tape<unsigned long> _t30 = {};
        clad::tape<int> _t31 = {};
        clad::array<double> _delta_Ap(_d_Ap.size());
        clad::array<double> _EERepl_Ap0(_d_Ap.size());
        for (int i0 = 0; i0 < _d_Ap.size(); i0++)
        {
            _EERepl_Ap0[i0] = Ap[i0];
        }
        clad::tape<double> _t33 = {};
        clad::tape<int> _t34 = {};
        clad::tape<int> _t36 = {};
        clad::tape<double> _t38 = {};
        clad::tape<int> _t39 = {};
        clad::tape<int> _t41 = {};
        clad::tape<int> _t43 = {};
        clad::tape<double> _EERepl_Ap1 = {};
        double _EERepl_alp3;
        double _EERepl_bta2;
        bool _cond4;
        unsigned long _t45;
        clad::tape<int> _t46 = {};
        clad::array<double> _delta_r(_d_r.size());
        clad::array<double> _EERepl_r0(_d_r.size());
        for (int i0 = 0; i0 < _d_r.size(); i0++)
        {
            _EERepl_r0[i0] = r[i0];
        }
        clad::tape<int> _t48 = {};
        clad::tape<double> _t50 = {};
        clad::tape<double> _t51 = {};
        clad::tape<int> _t52 = {};
        clad::tape<double> _EERepl_r1 = {};
        bool _cond5;
        unsigned long _t54;
        clad::tape<int> _t55 = {};
        clad::tape<double> _t57 = {};
        clad::tape<double> _t58 = {};
        clad::tape<int> _t59 = {};
        clad::tape<int> _t61 = {};
        clad::tape<double> _EERepl_r2 = {};
        unsigned long _t63;
        clad::tape<int> _t64 = {};
        clad::tape<double> _t66 = {};
        clad::tape<double> _t67 = {};
        clad::tape<int> _t68 = {};
        clad::tape<double> _t70 = {};
        clad::tape<double> _t71 = {};
        clad::tape<int> _t72 = {};
        clad::tape<double> _EERepl_r3 = {};
        double _EERepl_rtrans1;
        unsigned long _t74;
        clad::tape<double> _t75 = {};
        clad::tape<int> _t76 = {};
        clad::tape<double> _t78 = {};
        clad::tape<int> _t79 = {};
        clad::tape<double> _EERepl_rtrans2 = {};
        double _EERepl_normr1 = normr;
        double _t81;
        unsigned long _t82;
        int _d_k = 0;
        clad::tape<bool> _t84 = {};
        clad::tape<double> _EERepl_alp4 = {};
        clad::tape<double> _EERepl_bta3 = {};
        clad::tape<bool> _t86 = {};
        clad::tape<unsigned long> _t87 = {};
        clad::tape<int> _t88 = {};
        clad::tape<int> _t90 = {};
        clad::tape<double> _t92 = {};
        clad::tape<double> _t93 = {};
        clad::tape<int> _t94 = {};
        clad::tape<double> _EERepl_p4 = {};
        clad::tape<bool> _t97 = {};
        clad::tape<unsigned long> _t98 = {};
        clad::tape<int> _t99 = {};
        clad::tape<double> _t101 = {};
        clad::tape<double> _t102 = {};
        clad::tape<int> _t103 = {};
        clad::tape<int> _t105 = {};
        clad::tape<double> _EERepl_p5 = {};
        clad::tape<unsigned long> _t107 = {};
        clad::tape<int> _t108 = {};
        clad::tape<double> _t110 = {};
        clad::tape<double> _t111 = {};
        clad::tape<int> _t112 = {};
        clad::tape<double> _t114 = {};
        clad::tape<double> _t115 = {};
        clad::tape<int> _t116 = {};
        clad::tape<double> _EERepl_p6 = {};
        clad::tape<double> _EERepl_rtrans3 = {};
        clad::tape<unsigned long> _t118 = {};
        clad::tape<double> _t119 = {};
        clad::tape<int> _t120 = {};
        clad::tape<double> _t122 = {};
        clad::tape<int> _t123 = {};
        clad::tape<double> _EERepl_rtrans4 = {};
        clad::tape<double> _t125 = {};
        clad::tape<double> _t126 = {};
        double _d_beta = 0;
        double _delta_beta = 0;
        clad::tape<double> _EERepl_beta0 = {};
        clad::tape<double> _EERepl_alp5 = {};
        clad::tape<bool> _t128 = {};
        clad::tape<unsigned long> _t129 = {};
        clad::tape<int> _t130 = {};
        clad::tape<int> _t132 = {};
        clad::tape<double> _t134 = {};
        clad::tape<double> _t135 = {};
        clad::tape<int> _t136 = {};
        clad::tape<double> _EERepl_p7 = {};
        clad::tape<bool> _t139 = {};
        clad::tape<unsigned long> _t140 = {};
        clad::tape<int> _t141 = {};
        clad::tape<double> _t143 = {};
        clad::tape<double> _t144 = {};
        clad::tape<int> _t145 = {};
        clad::tape<int> _t147 = {};
        clad::tape<double> _EERepl_p8 = {};
        clad::tape<unsigned long> _t149 = {};
        clad::tape<int> _t150 = {};
        clad::tape<double> _t152 = {};
        clad::tape<double> _t153 = {};
        clad::tape<int> _t154 = {};
        clad::tape<double> _t156 = {};
        clad::tape<double> _t157 = {};
        clad::tape<int> _t158 = {};
        clad::tape<double> _EERepl_p9 = {};
        double _EERepl_normr2 = normr;
        clad::tape<double> _t160 = {};
        clad::tape<unsigned long> _t161 = {};
        clad::tape<unsigned long> _t162 = {};
        clad::tape<int> _t163 = {};
        clad::tape<double> _t165 = {};
        clad::tape<int> _t166 = {};
        clad::tape<int> _t168 = {};
        clad::tape<double> _t170 = {};
        clad::tape<int> _t171 = {};
        clad::tape<int> _t173 = {};
        clad::tape<int> _t175 = {};
        clad::tape<double> _EERepl_Ap2 = {};
        double _d_alpha = 0;
        double _delta_alpha = 0;
        clad::tape<double> _EERepl_alpha0 = {};
        clad::tape<double> _EERepl_alpha1 = {};
        clad::tape<bool> _t178 = {};
        clad::tape<unsigned long> _t179 = {};
        clad::tape<double> _t180 = {};
        clad::tape<int> _t181 = {};
        clad::tape<double> _t183 = {};
        clad::tape<int> _t184 = {};
        clad::tape<double> _EERepl_alpha2 = {};
        clad::tape<unsigned long> _t186 = {};
        clad::tape<double> _t187 = {};
        clad::tape<int> _t188 = {};
        clad::tape<double> _t190 = {};
        clad::tape<int> _t191 = {};
        clad::tape<double> _EERepl_alpha3 = {};
        clad::tape<double> _t193 = {};
        clad::tape<double> _t194 = {};
        clad::tape<double> _EERepl_alpha4 = {};
        clad::tape<double> _EERepl_alp6 = {};
        clad::tape<bool> _t196 = {};
        clad::tape<unsigned long> _t197 = {};
        clad::tape<int> _t198 = {};
        clad::array<double> _delta_x(_d_x.size());
        clad::array<double> _EERepl_x0(_d_x.size());
        for (int i0 = 0; i0 < _d_x.size(); i0++)
        {
            _EERepl_x0[i0] = x[i0];
        }
        clad::tape<int> _t200 = {};
        clad::tape<double> _t202 = {};
        clad::tape<double> _t203 = {};
        clad::tape<int> _t204 = {};
        clad::tape<double> _EERepl_x1 = {};
        clad::tape<bool> _t207 = {};
        clad::tape<unsigned long> _t208 = {};
        clad::tape<int> _t209 = {};
        clad::tape<double> _t211 = {};
        clad::tape<double> _t212 = {};
        clad::tape<int> _t213 = {};
        clad::tape<int> _t215 = {};
        clad::tape<double> _EERepl_x2 = {};
        clad::tape<unsigned long> _t217 = {};
        clad::tape<int> _t218 = {};
        clad::tape<double> _t220 = {};
        clad::tape<double> _t221 = {};
        clad::tape<int> _t222 = {};
        clad::tape<double> _t224 = {};
        clad::tape<double> _t225 = {};
        clad::tape<int> _t226 = {};
        clad::tape<double> _EERepl_x3 = {};
        clad::tape<double> _EERepl_alp7 = {};
        clad::tape<double> _EERepl_bta4 = {};
        clad::tape<bool> _t229 = {};
        clad::tape<unsigned long> _t230 = {};
        clad::tape<int> _t231 = {};
        clad::tape<int> _t233 = {};
        clad::tape<double> _t235 = {};
        clad::tape<double> _t236 = {};
        clad::tape<int> _t237 = {};
        clad::tape<double> _EERepl_r4 = {};
        clad::tape<bool> _t240 = {};
        clad::tape<unsigned long> _t241 = {};
        clad::tape<int> _t242 = {};
        clad::tape<double> _t244 = {};
        clad::tape<double> _t245 = {};
        clad::tape<int> _t246 = {};
        clad::tape<int> _t248 = {};
        clad::tape<double> _EERepl_r5 = {};
        clad::tape<unsigned long> _t250 = {};
        clad::tape<int> _t251 = {};
        clad::tape<double> _t253 = {};
        clad::tape<double> _t254 = {};
        clad::tape<int> _t255 = {};
        clad::tape<double> _t257 = {};
        clad::tape<double> _t258 = {};
        clad::tape<int> _t259 = {};
        clad::tape<double> _EERepl_r6 = {};
        double _d_residual = 0;
        double _delta_residual = 0;
        double _EERepl_residual0;
        unsigned long _t261;
        clad::tape<int> _t262 = {};
        clad::tape<int> _t264 = {};
        clad::tape<double> _t266 = {};
        double _d_diff = 0;
        double _delta_diff = 0;
        clad::tape<double> _EERepl_diff0 = {};
        clad::tape<bool> _t269 = {};
        int nrow = A.local_nrow;
        int ncol = A.local_ncol;
        normr = 0.;
        double rtrans = 0.;
        _EERepl_rtrans0 = rtrans;
        double oldrtrans = 0.;
        _EERepl_oldrtrans0 = oldrtrans;
        double alp, bta;
        _EERepl_alp1 = alp;
        _EERepl_alp0 = alp;
        int rank = 0;
        int print_freq = max_iter / 10;
        _cond0 = print_freq > 50;
        if (_cond0)
            print_freq = 50;
        _cond1 = print_freq < 1;
        if (_cond1)
            print_freq = 1;
        alp = 1.;
        _EERepl_alp2 = alp;
        bta = 0.;
        _EERepl_bta1 = bta;
        _cond2 = alp == 1.;
        if (_cond2)
        {
            _t0 = 0;
            for (int i = 0; i < nrow; i++)
            {
                _t0++;
                p[clad::push(_t1, i)] = x[clad::push(_t3, i)] + clad::push(_t6, bta) * clad::push(_t5, x[clad::push(_t7, i)]);
                clad::push(_EERepl_p1, p[clad::push(_t1, i)]);
            }
        }
        else
        {
            _cond3 = bta == 1.;
            if (_cond3)
            {
                _t9 = 0;
                for (int i = 0; i < nrow; i++)
                {
                    _t9++;
                    p[clad::push(_t10, i)] = clad::push(_t13, alp) * clad::push(_t12, x[clad::push(_t14, i)]) + x[clad::push(_t16, i)];
                    clad::push(_EERepl_p2, p[clad::push(_t10, i)]);
                }
            }
            else
            {
                _t18 = 0;
                for (int i = 0; i < nrow; i++)
                {
                    _t18++;
                    p[clad::push(_t19, i)] = clad::push(_t22, alp) * clad::push(_t21, x[clad::push(_t23, i)]) + clad::push(_t26, bta) * clad::push(_t25, x[clad::push(_t27, i)]);
                    clad::push(_EERepl_p3, p[clad::push(_t19, i)]);
                }
            }
        }
        _t29 = 0;
        for (int i = 0; i < nrow; i++)
        {
            _t29++;
            clad::push(_t30, 0UL);
            for (int j = 0; j < A.nnz_in_row[i]; j++)
            {
                clad::back(_t30)++;
                Ap[clad::push(_t31, i)] += clad::push(_t38, A.ptr_to_vals_in_row[clad::push(_t34, i)][clad::push(_t36, j)]) * clad::push(_t33, p[clad::push(_t43, A.ptr_to_inds_in_row[clad::push(_t39, i)][clad::push(_t41, j)])]);
                clad::push(_EERepl_Ap1, Ap[clad::push(_t31, i)]);
            }
        }
        alp = 1.;
        _EERepl_alp3 = alp;
        bta = -1.;
        _EERepl_bta2 = bta;
        _cond4 = alp == 1.;
        if (_cond4)
        {
            _t45 = 0;
            for (int i = 0; i < nrow; i++)
            {
                _t45++;
                r[clad::push(_t46, i)] = b[clad::push(_t48, i)] + clad::push(_t51, bta) * clad::push(_t50, Ap[clad::push(_t52, i)]);
                clad::push(_EERepl_r1, r[clad::push(_t46, i)]);
            }
        }
        else
        {
            _cond5 = bta == 1.;
            if (_cond5)
            {
                _t54 = 0;
                for (int i = 0; i < nrow; i++)
                {
                    _t54++;
                    r[clad::push(_t55, i)] = clad::push(_t58, alp) * clad::push(_t57, b[clad::push(_t59, i)]) + Ap[clad::push(_t61, i)];
                    clad::push(_EERepl_r2, r[clad::push(_t55, i)]);
                }
            }
            else
            {
                _t63 = 0;
                for (int i = 0; i < nrow; i++)
                {
                    _t63++;
                    r[clad::push(_t64, i)] = clad::push(_t67, alp) * clad::push(_t66, b[clad::push(_t68, i)]) + clad::push(_t71, bta) * clad::push(_t70, Ap[clad::push(_t72, i)]);
                    clad::push(_EERepl_r3, r[clad::push(_t64, i)]);
                }
            }
        }
        rtrans = 0.;
        _EERepl_rtrans1 = rtrans;
        _t74 = 0;
        for (int i = 0; i < nrow; i++)
        {
            _t74++;
            rtrans += clad::push(_t78, r[clad::push(_t76, i)]) * clad::push(_t75, r[clad::push(_t79, i)]);
            clad::push(_EERepl_rtrans2, rtrans);
        }
        _t81 = rtrans;
        normr = sqrt(_t81);
        _t82 = 0;
        for (int k = 1; k < max_iter && normr > tolerance; k++)
        {
            _t82++;
            bool _t83 = k == 1;
            {
                if (_t83)
                {
                    alp = 1.;
                    clad::push(_EERepl_alp4, alp);
                    bta = 0.;
                    clad::push(_EERepl_bta3, bta);
                    bool _t85 = alp == 1.;
                    {
                        if (_t85)
                        {
                            clad::push(_t87, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t87)++;
                                p[clad::push(_t88, i)] = r[clad::push(_t90, i)] + clad::push(_t93, bta) * clad::push(_t92, r[clad::push(_t94, i)]);
                                clad::push(_EERepl_p4, p[clad::push(_t88, i)]);
                            }
                        }
                        else
                        {
                            bool _t96 = bta == 1.;
                            {
                                if (_t96)
                                {
                                    clad::push(_t98, 0UL);
                                    for (int i = 0; i < nrow; i++)
                                    {
                                        clad::back(_t98)++;
                                        p[clad::push(_t99, i)] = clad::push(_t102, alp) * clad::push(_t101, r[clad::push(_t103, i)]) + r[clad::push(_t105, i)];
                                        clad::push(_EERepl_p5, p[clad::push(_t99, i)]);
                                    }
                                }
                                else
                                {
                                    clad::push(_t107, 0UL);
                                    for (int i = 0; i < nrow; i++)
                                    {
                                        clad::back(_t107)++;
                                        p[clad::push(_t108, i)] = clad::push(_t111, alp) * clad::push(_t110, r[clad::push(_t112, i)]) + clad::push(_t115, bta) * clad::push(_t114, r[clad::push(_t116, i)]);
                                        clad::push(_EERepl_p6, p[clad::push(_t108, i)]);
                                    }
                                }
                                clad::push(_t97, _t96);
                            }
                        }
                        clad::push(_t86, _t85);
                    }
                }
                else
                {
                    oldrtrans = rtrans;
                    rtrans = 0.;
                    clad::push(_EERepl_rtrans3, rtrans);
                    clad::push(_t118, 0UL);
                    for (int i = 0; i < nrow; i++)
                    {
                        clad::back(_t118)++;
                        rtrans += clad::push(_t122, r[clad::push(_t120, i)]) * clad::push(_t119, r[clad::push(_t123, i)]);
                        clad::push(_EERepl_rtrans4, rtrans);
                    }
                    double beta = clad::push(_t126, rtrans) / clad::push(_t125, oldrtrans);
                    clad::push(_EERepl_beta0, beta);
                    alp = 1.;
                    clad::push(_EERepl_alp5, alp);
                    bta = beta;
                    bool _t127 = alp == 1.;
                    {
                        if (_t127)
                        {
                            clad::push(_t129, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t129)++;
                                p[clad::push(_t130, i)] = r[clad::push(_t132, i)] + clad::push(_t135, bta) * clad::push(_t134, p[clad::push(_t136, i)]);
                                clad::push(_EERepl_p7, p[clad::push(_t130, i)]);
                            }
                        }
                        else
                        {
                            bool _t138 = bta == 1.;
                            {
                                if (_t138)
                                {
                                    clad::push(_t140, 0UL);
                                    for (int i = 0; i < nrow; i++)
                                    {
                                        clad::back(_t140)++;
                                        p[clad::push(_t141, i)] = clad::push(_t144, alp) * clad::push(_t143, r[clad::push(_t145, i)]) + p[clad::push(_t147, i)];
                                        clad::push(_EERepl_p8, p[clad::push(_t141, i)]);
                                    }
                                }
                                else
                                {
                                    clad::push(_t149, 0UL);
                                    for (int i = 0; i < nrow; i++)
                                    {
                                        clad::back(_t149)++;
                                        p[clad::push(_t150, i)] = clad::push(_t153, alp) * clad::push(_t152, r[clad::push(_t154, i)]) + clad::push(_t157, bta) * clad::push(_t156, p[clad::push(_t158, i)]);
                                        clad::push(_EERepl_p9, p[clad::push(_t150, i)]);
                                    }
                                }
                                clad::push(_t139, _t138);
                            }
                        }
                        clad::push(_t128, _t127);
                    }
                }
                clad::push(_t84, _t83);
            }
            normr = sqrt(clad::push(_t160, rtrans));
            clad::push(_t161, 0UL);
            for (int i = 0; i < nrow; i++)
            {
                clad::back(_t161)++;
                clad::push(_t162, 0UL);
                for (int j = 0; j < A.nnz_in_row[i]; j++)
                {
                    clad::back(_t162)++;
                    Ap[clad::push(_t163, i)] += clad::push(_t170, A.ptr_to_vals_in_row[clad::push(_t166, i)][clad::push(_t168, j)]) * clad::push(_t165, p[clad::push(_t175, A.ptr_to_inds_in_row[clad::push(_t171, i)][clad::push(_t173, j)])]);
                    clad::push(_EERepl_Ap2, Ap[clad::push(_t163, i)]);
                }
            }
            double alpha = 0.;
            clad::push(_EERepl_alpha0, alpha);
            alpha = 0.;
            clad::push(_EERepl_alpha1, alpha);
            bool _t177 = Ap == p;
            {
                if (_t177)
                {
                    clad::push(_t179, 0UL);
                    for (int i = 0; i < nrow; i++)
                    {
                        clad::back(_t179)++;
                        alpha += clad::push(_t183, p[clad::push(_t181, i)]) * clad::push(_t180, p[clad::push(_t184, i)]);
                        clad::push(_EERepl_alpha2, alpha);
                    }
                }
                else
                {
                    clad::push(_t186, 0UL);
                    for (int i = 0; i < nrow; i++)
                    {
                        clad::back(_t186)++;
                        alpha += clad::push(_t190, p[clad::push(_t188, i)]) * clad::push(_t187, Ap[clad::push(_t191, i)]);
                        clad::push(_EERepl_alpha3, alpha);
                    }
                }
                clad::push(_t178, _t177);
            }
            alpha = clad::push(_t194, rtrans) / clad::push(_t193, alpha);
            clad::push(_EERepl_alpha4, alpha);
            alp = 1.;
            clad::push(_EERepl_alp6, alp);
            bta = alpha;
            bool _t195 = alp == 1.;
            {
                if (_t195)
                {
                    clad::push(_t197, 0UL);
                    for (int i = 0; i < nrow; i++)
                    {
                        clad::back(_t197)++;
                        x[clad::push(_t198, i)] = x[clad::push(_t200, i)] + clad::push(_t203, bta) * clad::push(_t202, p[clad::push(_t204, i)]);
                        clad::push(_EERepl_x1, x[clad::push(_t198, i)]);
                    }
                }
                else
                {
                    bool _t206 = bta == 1.;
                    {
                        if (_t206)
                        {
                            clad::push(_t208, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t208)++;
                                x[clad::push(_t209, i)] = clad::push(_t212, alp) * clad::push(_t211, x[clad::push(_t213, i)]) + p[clad::push(_t215, i)];
                                clad::push(_EERepl_x2, x[clad::push(_t209, i)]);
                            }
                        }
                        else
                        {
                            clad::push(_t217, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t217)++;
                                x[clad::push(_t218, i)] = clad::push(_t221, alp) * clad::push(_t220, x[clad::push(_t222, i)]) + clad::push(_t225, bta) * clad::push(_t224, p[clad::push(_t226, i)]);
                                clad::push(_EERepl_x3, x[clad::push(_t218, i)]);
                            }
                        }
                        clad::push(_t207, _t206);
                    }
                }
                clad::push(_t196, _t195);
            }
            alp = 1.;
            clad::push(_EERepl_alp7, alp);
            bta = -alpha;
            clad::push(_EERepl_bta4, bta);
            bool _t228 = alp == 1.;
            {
                if (_t228)
                {
                    clad::push(_t230, 0UL);
                    for (int i = 0; i < nrow; i++)
                    {
                        clad::back(_t230)++;
                        r[clad::push(_t231, i)] = r[clad::push(_t233, i)] + clad::push(_t236, bta) * clad::push(_t235, Ap[clad::push(_t237, i)]);
                        clad::push(_EERepl_r4, r[clad::push(_t231, i)]);
                    }
                }
                else
                {
                    bool _t239 = bta == 1.;
                    {
                        if (_t239)
                        {
                            clad::push(_t241, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t241)++;
                                r[clad::push(_t242, i)] = clad::push(_t245, alp) * clad::push(_t244, r[clad::push(_t246, i)]) + Ap[clad::push(_t248, i)];
                                clad::push(_EERepl_r5, r[clad::push(_t242, i)]);
                            }
                        }
                        else
                        {
                            clad::push(_t250, 0UL);
                            for (int i = 0; i < nrow; i++)
                            {
                                clad::back(_t250)++;
                                r[clad::push(_t251, i)] = clad::push(_t254, alp) * clad::push(_t253, r[clad::push(_t255, i)]) + clad::push(_t258, bta) * clad::push(_t257, Ap[clad::push(_t259, i)]);
                                clad::push(_EERepl_r6, r[clad::push(_t251, i)]);
                            }
                        }
                        clad::push(_t240, _t239);
                    }
                }
                clad::push(_t229, _t228);
            }
            niters = k;
        }
        double residual = 0;
        _EERepl_residual0 = residual;
        _t261 = 0;
        for (int i = 0; i < nrow; i++)
        {
            _t261++;
            double diff = fabs(clad::push(_t266, x[clad::push(_t262, i)] - xexact[clad::push(_t264, i)]));
            clad::push(_EERepl_diff0, diff);
            bool _t268 = diff > residual;
            {
                if (_t268)
                    residual = diff;
                clad::push(_t269, _t268);
            }
        }
        double HPCCG_return = residual;
        goto _label0;
    _label0:
        _d_residual += 1;
        for (; _t261; _t261--)
        {
            if (clad::pop(_t269))
            {
                double _r_d47 = _d_residual;
                _d_diff += _r_d47;
                _d_residual -= _r_d47;
            }
            {
                double _r124 = _d_diff;
                int _t263 = clad::pop(_t262);
                _d_x[_t263] += _r124;
                int _t265 = clad::pop(_t264);
                _d_xexact[_t265] += -_r124;
                _d_diff = 0;
                double _r125 = clad::pop(_EERepl_diff0);
                _delta_diff += clad::getErrorVal(_d_diff, _r125, "diff");;
            }
        }
        _delta_residual += clad::getErrorVal(_d_residual, _EERepl_residual0, "residual");
        for (; _t82; _t82--)
        {
            {
                int _r_d46 = *_d_niters;
                _d_k += _r_d46;
                *_d_niters -= _r_d46;
                *_d_niters;
            }
            if (clad::pop(_t229))
            {
                {
                    for (; clad::back(_t230); clad::back(_t230)--)
                    {
                        int _t232 = clad::pop(_t231);
                        double _r_d43 = _d_r[_t232];
                        int _t234 = clad::pop(_t233);
                        _d_r[_t234] += _r_d43;
                        double _r110 = _r_d43 * clad::pop(_t235);
                        _d_bta += _r110;
                        double _r111 = clad::pop(_t236) * _r_d43;
                        int _t238 = clad::pop(_t237);
                        _d_Ap[_t238] += _r111;
                        double _r112 = clad::pop(_EERepl_r4);
                        int _r113 = clad::pop(_t231);
                        _delta_r[_r113] += clad::getErrorVal(_r_d43, _r112, "r");
                        _final_error += _delta_r[_r113];
                        _d_r[_t232] -= _r_d43;
                        _d_r[_t232];
                    }
                    clad::pop(_t230);
                }
            }
            else if (clad::pop(_t240))
            {
                {
                    for (; clad::back(_t241); clad::back(_t241)--)
                    {
                        int _t243 = clad::pop(_t242);
                        double _r_d44 = _d_r[_t243];
                        double _r114 = _r_d44 * clad::pop(_t244);
                        _d_alp += _r114;
                        double _r115 = clad::pop(_t245) * _r_d44;
                        int _t247 = clad::pop(_t246);
                        _d_r[_t247] += _r115;
                        int _t249 = clad::pop(_t248);
                        _d_Ap[_t249] += _r_d44;
                        double _r116 = clad::pop(_EERepl_r5);
                        int _r117 = clad::pop(_t242);
                        _delta_r[_r117] += clad::getErrorVal(_r_d44, _r116, "r");
                        _final_error += _delta_r[_r117];
                        _d_r[_t243] -= _r_d44;
                        _d_r[_t243];
                    }
                    clad::pop(_t241);
                }
            }
            else
            {
                {
                    for (; clad::back(_t250); clad::back(_t250)--)
                    {
                        int _t252 = clad::pop(_t251);
                        double _r_d45 = _d_r[_t252];
                        double _r118 = _r_d45 * clad::pop(_t253);
                        _d_alp += _r118;
                        double _r119 = clad::pop(_t254) * _r_d45;
                        int _t256 = clad::pop(_t255);
                        _d_r[_t256] += _r119;
                        double _r120 = _r_d45 * clad::pop(_t257);
                        _d_bta += _r120;
                        double _r121 = clad::pop(_t258) * _r_d45;
                        int _t260 = clad::pop(_t259);
                        _d_Ap[_t260] += _r121;
                        double _r122 = clad::pop(_EERepl_r6);
                        int _r123 = clad::pop(_t251);
                        _delta_r[_r123] += clad::getErrorVal(_r_d45, _r122, "r");
                        _final_error += _delta_r[_r123];
                        _d_r[_t252] -= _r_d45;
                        _d_r[_t252];
                    }
                    clad::pop(_t250);
                }
            }
            {
                double _r_d42 = _d_bta;
                _d_alpha += -_r_d42;
                double _r109 = clad::pop(_EERepl_bta4);
                _delta_bta += clad::getErrorVal(_r_d42, _r109, "bta");
                _d_bta -= _r_d42;
            }
            {
                double _r_d41 = _d_alp;
                double _r108 = clad::pop(_EERepl_alp7);
                _delta_alp += clad::getErrorVal(_r_d41, _r108, "alp");
                _d_alp -= _r_d41;
            }
            if (clad::pop(_t196))
            {
                {
                    for (; clad::back(_t197); clad::back(_t197)--)
                    {
                        int _t199 = clad::pop(_t198);
                        double _r_d38 = _d_x[_t199];
                        int _t201 = clad::pop(_t200);
                        _d_x[_t201] += _r_d38;
                        double _r94 = _r_d38 * clad::pop(_t202);
                        _d_bta += _r94;
                        double _r95 = clad::pop(_t203) * _r_d38;
                        int _t205 = clad::pop(_t204);
                        _d_p[_t205] += _r95;
                        double _r96 = clad::pop(_EERepl_x1);
                        int _r97 = clad::pop(_t198);
                        _delta_x[_r97] += clad::getErrorVal(_r_d38, _r96, "x");
                        _final_error += _delta_x[_r97];
                        _d_x[_t199] -= _r_d38;
                        _d_x[_t199];
                    }
                    clad::pop(_t197);
                }
            }
            else if (clad::pop(_t207))
            {
                {
                    for (; clad::back(_t208); clad::back(_t208)--)
                    {
                        int _t210 = clad::pop(_t209);
                        double _r_d39 = _d_x[_t210];
                        double _r98 = _r_d39 * clad::pop(_t211);
                        _d_alp += _r98;
                        double _r99 = clad::pop(_t212) * _r_d39;
                        int _t214 = clad::pop(_t213);
                        _d_x[_t214] += _r99;
                        int _t216 = clad::pop(_t215);
                        _d_p[_t216] += _r_d39;
                        double _r100 = clad::pop(_EERepl_x2);
                        int _r101 = clad::pop(_t209);
                        _delta_x[_r101] += clad::getErrorVal(_r_d39, _r100, "x");
                        _final_error += _delta_x[_r101];
                        _d_x[_t210] -= _r_d39;
                        _d_x[_t210];
                    }
                    clad::pop(_t208);
                }
            }
            else
            {
                {
                    for (; clad::back(_t217); clad::back(_t217)--)
                    {
                        int _t219 = clad::pop(_t218);
                        double _r_d40 = _d_x[_t219];
                        double _r102 = _r_d40 * clad::pop(_t220);
                        _d_alp += _r102;
                        double _r103 = clad::pop(_t221) * _r_d40;
                        int _t223 = clad::pop(_t222);
                        _d_x[_t223] += _r103;
                        double _r104 = _r_d40 * clad::pop(_t224);
                        _d_bta += _r104;
                        double _r105 = clad::pop(_t225) * _r_d40;
                        int _t227 = clad::pop(_t226);
                        _d_p[_t227] += _r105;
                        double _r106 = clad::pop(_EERepl_x3);
                        int _r107 = clad::pop(_t218);
                        _delta_x[_r107] += clad::getErrorVal(_r_d40, _r106, "x");
                        _final_error += _delta_x[_r107];
                        _d_x[_t219] -= _r_d40;
                        _d_x[_t219];
                    }
                    clad::pop(_t217);
                }
            }
            {
                double _r_d37 = _d_bta;
                _d_alpha += _r_d37;
                _d_bta -= _r_d37;
            }
            {
                double _r_d36 = _d_alp;
                double _r93 = clad::pop(_EERepl_alp6);
                _delta_alp += clad::getErrorVal(_r_d36, _r93, "alp");
                _d_alp -= _r_d36;
            }
            {
                double _r_d35 = _d_alpha;
                double _r89 = clad::pop(_t193);
                double _r90 = _r_d35 / _r89;
                _d_rtrans += _r90;
                double _r91 = _r_d35 * -clad::pop(_t194) / (_r89 * _r89);
                _d_alpha += _r91;
                double _r92 = clad::pop(_EERepl_alpha4);
                _delta_alpha += clad::getErrorVal(_r_d35, _r92, "alpha");
                _d_alpha -= _r_d35;
            }
            if (clad::pop(_t178))
            {
                for (; clad::back(_t179); clad::back(_t179)--)
                {
                    double _r_d33 = _d_alpha;
                    _d_alpha += _r_d33;
                    double _r83 = _r_d33 * clad::pop(_t180);
                    int _t182 = clad::pop(_t181);
                    _d_p[_t182] += _r83;
                    double _r84 = clad::pop(_t183) * _r_d33;
                    int _t185 = clad::pop(_t184);
                    _d_p[_t185] += _r84;
                    double _r85 = clad::pop(_EERepl_alpha2);
                    _delta_alpha += clad::getErrorVal(_r_d33, _r85, "alpha");
                    _d_alpha -= _r_d33;
                }
                clad::pop(_t179);
            }
            else
            {
                for (; clad::back(_t186); clad::back(_t186)--)
                {
                    double _r_d34 = _d_alpha;
                    _d_alpha += _r_d34;
                    double _r86 = _r_d34 * clad::pop(_t187);
                    int _t189 = clad::pop(_t188);
                    _d_p[_t189] += _r86;
                    double _r87 = clad::pop(_t190) * _r_d34;
                    int _t192 = clad::pop(_t191);
                    _d_Ap[_t192] += _r87;
                    double _r88 = clad::pop(_EERepl_alpha3);
                    _delta_alpha += clad::getErrorVal(_r_d34, _r88, "alpha");
                    _d_alpha -= _r_d34;
                }
                clad::pop(_t186);
            }
            {
                int _t174 = clad::pop(_t173);
                int _t172 = clad::pop(_t171);
                double _r_d32 = _d_alpha;
                double _r82 = clad::pop(_EERepl_alpha1);
                _delta_alpha += clad::getErrorVal(_r_d32, _r82, "alpha");
                _d_alpha -= _r_d32;
            }
            {
                _d_alpha = 0;
                double _r81 = clad::pop(_EERepl_alpha0);
                _delta_alpha += clad::getErrorVal(_d_alpha, _r81, "alpha");
            }
            {
                for (; clad::back(_t161); clad::back(_t161)--)
                {
                    {
                        for (; clad::back(_t162); clad::back(_t162)--)
                        {
                            int _t164 = clad::pop(_t163);
                            double _r_d31 = _d_Ap[_t164];
                            _d_Ap[_t164] += _r_d31;
                            double _r77 = _r_d31 * clad::pop(_t165);
                            int _t167 = clad::pop(_t166);
                            int _t169 = clad::pop(_t168);
                            double _r78 = clad::pop(_t170) * _r_d31;
                            int _t176 = clad::pop(_t175);
                            _d_p[_t176] += _r78;
                            double _r79 = clad::pop(_EERepl_Ap2);
                            int _r80 = clad::pop(_t163);
                            _delta_Ap[_r80] += clad::getErrorVal(_r_d31, _r79, "Ap");
                            _final_error += _delta_Ap[_r80];
                            _d_Ap[_t164] -= _r_d31;
                            _d_Ap[_t164];
                        }
                        clad::pop(_t162);
                    }
                }
                clad::pop(_t161);
            }
            {
                double _r_d30 = *_d_normr;
                double _r76 = _r_d30 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t160), 1.).pushforward;
                _d_rtrans += _r76;
                *_d_normr -= _r_d30;
                *_d_normr;
            }
            if (clad::pop(_t84))
            {
                if (clad::pop(_t86))
                {
                    {
                        for (; clad::back(_t87); clad::back(_t87)--)
                        {
                            int _t89 = clad::pop(_t88);
                            double _r_d19 = _d_p[_t89];
                            int _t91 = clad::pop(_t90);
                            _d_r[_t91] += _r_d19;
                            double _r39 = _r_d19 * clad::pop(_t92);
                            _d_bta += _r39;
                            double _r40 = clad::pop(_t93) * _r_d19;
                            int _t95 = clad::pop(_t94);
                            _d_r[_t95] += _r40;
                            double _r41 = clad::pop(_EERepl_p4);
                            int _r42 = clad::pop(_t88);
                            _delta_p[_r42] += clad::getErrorVal(_r_d19, _r41, "p");
                            _final_error += _delta_p[_r42];
                            _d_p[_t89] -= _r_d19;
                            _d_p[_t89];
                        }
                        clad::pop(_t87);
                    }
                }
                else if (clad::pop(_t97))
                {
                    {
                        for (; clad::back(_t98); clad::back(_t98)--)
                        {
                            int _t100 = clad::pop(_t99);
                            double _r_d20 = _d_p[_t100];
                            double _r43 = _r_d20 * clad::pop(_t101);
                            _d_alp += _r43;
                            double _r44 = clad::pop(_t102) * _r_d20;
                            int _t104 = clad::pop(_t103);
                            _d_r[_t104] += _r44;
                            int _t106 = clad::pop(_t105);
                            _d_r[_t106] += _r_d20;
                            double _r45 = clad::pop(_EERepl_p5);
                            int _r46 = clad::pop(_t99);
                            _delta_p[_r46] += clad::getErrorVal(_r_d20, _r45, "p");
                            _final_error += _delta_p[_r46];
                            _d_p[_t100] -= _r_d20;
                            _d_p[_t100];
                        }
                        clad::pop(_t98);
                    }
                }
                else
                {
                    {
                        for (; clad::back(_t107); clad::back(_t107)--)
                        {
                            int _t109 = clad::pop(_t108);
                            double _r_d21 = _d_p[_t109];
                            double _r47 = _r_d21 * clad::pop(_t110);
                            _d_alp += _r47;
                            double _r48 = clad::pop(_t111) * _r_d21;
                            int _t113 = clad::pop(_t112);
                            _d_r[_t113] += _r48;
                            double _r49 = _r_d21 * clad::pop(_t114);
                            _d_bta += _r49;
                            double _r50 = clad::pop(_t115) * _r_d21;
                            int _t117 = clad::pop(_t116);
                            _d_r[_t117] += _r50;
                            double _r51 = clad::pop(_EERepl_p6);
                            int _r52 = clad::pop(_t108);
                            _delta_p[_r52] += clad::getErrorVal(_r_d21, _r51, "p");
                            _final_error += _delta_p[_r52];
                            _d_p[_t109] -= _r_d21;
                            _d_p[_t109];
                        }
                        clad::pop(_t107);
                    }
                }
                {
                    double _r_d18 = _d_bta;
                    double _r38 = clad::pop(_EERepl_bta3);
                    _delta_bta += clad::getErrorVal(_r_d18, _r38, "bta");
                    _d_bta -= _r_d18;
                }
                {
                    double _r_d17 = _d_alp;
                    double _r37 = clad::pop(_EERepl_alp4);
                    _delta_alp += clad::getErrorVal(_r_d17, _r37, "alp");
                    _d_alp -= _r_d17;
                }
            }
            else
            {
                if (clad::pop(_t128))
                {
                    {
                        for (; clad::back(_t129); clad::back(_t129)--)
                        {
                            int _t131 = clad::pop(_t130);
                            double _r_d27 = _d_p[_t131];
                            int _t133 = clad::pop(_t132);
                            _d_r[_t133] += _r_d27;
                            double _r62 = _r_d27 * clad::pop(_t134);
                            _d_bta += _r62;
                            double _r63 = clad::pop(_t135) * _r_d27;
                            int _t137 = clad::pop(_t136);
                            _d_p[_t137] += _r63;
                            double _r64 = clad::pop(_EERepl_p7);
                            int _r65 = clad::pop(_t130);
                            _delta_p[_r65] += clad::getErrorVal(_r_d27, _r64, "p");
                            _final_error += _delta_p[_r65];
                            _d_p[_t131] -= _r_d27;
                            _d_p[_t131];
                        }
                        clad::pop(_t129);
                    }
                }
                else if (clad::pop(_t139))
                {
                    {
                        for (; clad::back(_t140); clad::back(_t140)--)
                        {
                            int _t142 = clad::pop(_t141);
                            double _r_d28 = _d_p[_t142];
                            double _r66 = _r_d28 * clad::pop(_t143);
                            _d_alp += _r66;
                            double _r67 = clad::pop(_t144) * _r_d28;
                            int _t146 = clad::pop(_t145);
                            _d_r[_t146] += _r67;
                            int _t148 = clad::pop(_t147);
                            _d_p[_t148] += _r_d28;
                            double _r68 = clad::pop(_EERepl_p8);
                            int _r69 = clad::pop(_t141);
                            _delta_p[_r69] += clad::getErrorVal(_r_d28, _r68, "p");
                            _final_error += _delta_p[_r69];
                            _d_p[_t142] -= _r_d28;
                            _d_p[_t142];
                        }
                        clad::pop(_t140);
                    }
                }
                else
                {
                    {
                        for (; clad::back(_t149); clad::back(_t149)--)
                        {
                            int _t151 = clad::pop(_t150);
                            double _r_d29 = _d_p[_t151];
                            double _r70 = _r_d29 * clad::pop(_t152);
                            _d_alp += _r70;
                            double _r71 = clad::pop(_t153) * _r_d29;
                            int _t155 = clad::pop(_t154);
                            _d_r[_t155] += _r71;
                            double _r72 = _r_d29 * clad::pop(_t156);
                            _d_bta += _r72;
                            double _r73 = clad::pop(_t157) * _r_d29;
                            int _t159 = clad::pop(_t158);
                            _d_p[_t159] += _r73;
                            double _r74 = clad::pop(_EERepl_p9);
                            int _r75 = clad::pop(_t150);
                            _delta_p[_r75] += clad::getErrorVal(_r_d29, _r74, "p");
                            _final_error += _delta_p[_r75];
                            _d_p[_t151] -= _r_d29;
                            _d_p[_t151];
                        }
                        clad::pop(_t149);
                    }
                }
                {
                    double _r_d26 = _d_bta;
                    _d_beta += _r_d26;
                    _d_bta -= _r_d26;
                }
                {
                    double _r_d25 = _d_alp;
                    double _r61 = clad::pop(_EERepl_alp5);
                    _delta_alp += clad::getErrorVal(_r_d25, _r61, "alp");
                    _d_alp -= _r_d25;
                }
                {
                    double _r57 = clad::pop(_t125);
                    double _r58 = _d_beta / _r57;
                    _d_rtrans += _r58;
                    double _r59 = _d_beta * -clad::pop(_t126) / (_r57 * _r57);
                    _d_oldrtrans += _r59;
                    _d_beta = 0;
                    double _r60 = clad::pop(_EERepl_beta0);
                    _delta_beta += clad::getErrorVal(_d_beta, _r60, "beta");
                }
                {
                    for (; clad::back(_t118); clad::back(_t118)--)
                    {
                        double _r_d24 = _d_rtrans;
                        _d_rtrans += _r_d24;
                        double _r54 = _r_d24 * clad::pop(_t119);
                        int _t121 = clad::pop(_t120);
                        _d_r[_t121] += _r54;
                        double _r55 = clad::pop(_t122) * _r_d24;
                        int _t124 = clad::pop(_t123);
                        _d_r[_t124] += _r55;
                        double _r56 = clad::pop(_EERepl_rtrans4);
                        _delta_rtrans += clad::getErrorVal(_r_d24, _r56, "rtrans");
                        _d_rtrans -= _r_d24;
                    }
                    clad::pop(_t118);
                }
                {
                    double _r_d23 = _d_rtrans;
                    double _r53 = clad::pop(_EERepl_rtrans3);
                    _delta_rtrans += clad::getErrorVal(_r_d23, _r53, "rtrans");
                    _d_rtrans -= _r_d23;
                }
                {
                    double _r_d22 = _d_oldrtrans;
                    _d_rtrans += _r_d22;
                    _d_oldrtrans -= _r_d22;
                }
            }
        }
        {
            double _r_d16 = *_d_normr;
            double _r36 = _r_d16 * clad::custom_derivatives::sqrt_pushforward(_t81, 1.).pushforward;
            _d_rtrans += _r36;
            *_d_normr -= _r_d16;
            *_d_normr;
        }
        for (; _t74; _t74--)
        {
            double _r_d15 = _d_rtrans;
            _d_rtrans += _r_d15;
            double _r33 = _r_d15 * clad::pop(_t75);
            int _t77 = clad::pop(_t76);
            _d_r[_t77] += _r33;
            double _r34 = clad::pop(_t78) * _r_d15;
            int _t80 = clad::pop(_t79);
            _d_r[_t80] += _r34;
            double _r35 = clad::pop(_EERepl_rtrans2);
            _delta_rtrans += clad::getErrorVal(_r_d15, _r35, "rtrans");
            _d_rtrans -= _r_d15;
        }
        {
            double _r_d14 = _d_rtrans;
            _delta_rtrans += clad::getErrorVal(_r_d14, _EERepl_rtrans1, "rtrans");
            _d_rtrans -= _r_d14;
        }
        if (_cond4)
        {
            for (; _t45; _t45--)
            {
                int _t47 = clad::pop(_t46);
                double _r_d11 = _d_r[_t47];
                int _t49 = clad::pop(_t48);
                _d_b[_t49] += _r_d11;
                double _r19 = _r_d11 * clad::pop(_t50);
                _d_bta += _r19;
                double _r20 = clad::pop(_t51) * _r_d11;
                int _t53 = clad::pop(_t52);
                _d_Ap[_t53] += _r20;
                double _r21 = clad::pop(_EERepl_r1);
                int _r22 = clad::pop(_t46);
                _delta_r[_r22] += clad::getErrorVal(_r_d11, _r21, "r");
                _final_error += _delta_r[_r22];
                _d_r[_t47] -= _r_d11;
                _d_r[_t47];
            }
        }
        else if (_cond5)
        {
            for (; _t54; _t54--)
            {
                int _t56 = clad::pop(_t55);
                double _r_d12 = _d_r[_t56];
                double _r23 = _r_d12 * clad::pop(_t57);
                _d_alp += _r23;
                double _r24 = clad::pop(_t58) * _r_d12;
                int _t60 = clad::pop(_t59);
                _d_b[_t60] += _r24;
                int _t62 = clad::pop(_t61);
                _d_Ap[_t62] += _r_d12;
                double _r25 = clad::pop(_EERepl_r2);
                int _r26 = clad::pop(_t55);
                _delta_r[_r26] += clad::getErrorVal(_r_d12, _r25, "r");
                _final_error += _delta_r[_r26];
                _d_r[_t56] -= _r_d12;
                _d_r[_t56];
            }
        }
        else
        {
            for (; _t63; _t63--)
            {
                int _t65 = clad::pop(_t64);
                double _r_d13 = _d_r[_t65];
                double _r27 = _r_d13 * clad::pop(_t66);
                _d_alp += _r27;
                double _r28 = clad::pop(_t67) * _r_d13;
                int _t69 = clad::pop(_t68);
                _d_b[_t69] += _r28;
                double _r29 = _r_d13 * clad::pop(_t70);
                _d_bta += _r29;
                double _r30 = clad::pop(_t71) * _r_d13;
                int _t73 = clad::pop(_t72);
                _d_Ap[_t73] += _r30;
                double _r31 = clad::pop(_EERepl_r3);
                int _r32 = clad::pop(_t64);
                _delta_r[_r32] += clad::getErrorVal(_r_d13, _r31, "r");
                _final_error += _delta_r[_r32];
                _d_r[_t65] -= _r_d13;
                _d_r[_t65];
            }
        }
        {
            double _r_d10 = _d_bta;
            _delta_bta += clad::getErrorVal(_r_d10, _EERepl_bta2, "bta");
            _d_bta -= _r_d10;
        }
        {
            int _t42 = clad::pop(_t41);
            int _t40 = clad::pop(_t39);
            double _r_d9 = _d_alp;
            _delta_alp += clad::getErrorVal(_r_d9, _EERepl_alp3, "alp");
            _d_alp -= _r_d9;
        }
        for (; _t29; _t29--)
        {
            {
                for (; clad::back(_t30); clad::back(_t30)--)
                {
                    int _t32 = clad::pop(_t31);
                    double _r_d8 = _d_Ap[_t32];
                    _d_Ap[_t32] += _r_d8;
                    double _r15 = _r_d8 * clad::pop(_t33);
                    int _t35 = clad::pop(_t34);
                    int _t37 = clad::pop(_t36);
                    double _r16 = clad::pop(_t38) * _r_d8;
                    int _t44 = clad::pop(_t43);
                    _d_p[_t44] += _r16;
                    double _r17 = clad::pop(_EERepl_Ap1);
                    int _r18 = clad::pop(_t31);
                    _delta_Ap[_r18] += clad::getErrorVal(_r_d8, _r17, "Ap");
                    _final_error += _delta_Ap[_r18];
                    _d_Ap[_t32] -= _r_d8;
                    _d_Ap[_t32];
                }
                clad::pop(_t30);
            }
        }
        if (_cond2)
        {
            for (; _t0; _t0--)
            {
                int _t2 = clad::pop(_t1);
                double _r_d5 = _d_p[_t2];
                int _t4 = clad::pop(_t3);
                _d_x[_t4] += _r_d5;
                double _r1 = _r_d5 * clad::pop(_t5);
                _d_bta += _r1;
                double _r2 = clad::pop(_t6) * _r_d5;
                int _t8 = clad::pop(_t7);
                _d_x[_t8] += _r2;
                double _r3 = clad::pop(_EERepl_p1);
                int _r4 = clad::pop(_t1);
                _delta_p[_r4] += clad::getErrorVal(_r_d5, _r3, "p");
                _final_error += _delta_p[_r4];
                _d_p[_t2] -= _r_d5;
                _d_p[_t2];
            }
        }
        else if (_cond3)
        {
            for (; _t9; _t9--)
            {
                int _t11 = clad::pop(_t10);
                double _r_d6 = _d_p[_t11];
                double _r5 = _r_d6 * clad::pop(_t12);
                _d_alp += _r5;
                double _r6 = clad::pop(_t13) * _r_d6;
                int _t15 = clad::pop(_t14);
                _d_x[_t15] += _r6;
                int _t17 = clad::pop(_t16);
                _d_x[_t17] += _r_d6;
                double _r7 = clad::pop(_EERepl_p2);
                int _r8 = clad::pop(_t10);
                _delta_p[_r8] += clad::getErrorVal(_r_d6, _r7, "p");
                _final_error += _delta_p[_r8];
                _d_p[_t11] -= _r_d6;
                _d_p[_t11];
            }
        }
        else
        {
            for (; _t18; _t18--)
            {
                int _t20 = clad::pop(_t19);
                double _r_d7 = _d_p[_t20];
                double _r9 = _r_d7 * clad::pop(_t21);
                _d_alp += _r9;
                double _r10 = clad::pop(_t22) * _r_d7;
                int _t24 = clad::pop(_t23);
                _d_x[_t24] += _r10;
                double _r11 = _r_d7 * clad::pop(_t25);
                _d_bta += _r11;
                double _r12 = clad::pop(_t26) * _r_d7;
                int _t28 = clad::pop(_t27);
                _d_x[_t28] += _r12;
                double _r13 = clad::pop(_EERepl_p3);
                int _r14 = clad::pop(_t19);
                _delta_p[_r14] += clad::getErrorVal(_r_d7, _r13, "p");
                _final_error += _delta_p[_r14];
                _d_p[_t20] -= _r_d7;
                _d_p[_t20];
            }
        }
        {
            double _r_d4 = _d_bta;
            _delta_bta += clad::getErrorVal(_r_d4, _EERepl_bta1, "bta");
            _d_bta -= _r_d4;
        }
        {
            double _r_d3 = _d_alp;
            _delta_alp += clad::getErrorVal(_r_d3, _EERepl_alp2, "alp");
            _d_alp -= _r_d3;
        }
        if (_cond1)
        {
            int _r_d2 = _d_print_freq;
            _d_print_freq -= _r_d2;
        }
        if (_cond0)
        {
            int _r_d1 = _d_print_freq;
            _d_print_freq -= _r_d1;
        }
        {
            int _r0 = _d_print_freq / 10;
            *_d_max_iter += _r0;
        }
        _delta_oldrtrans += clad::getErrorVal(_d_oldrtrans, _EERepl_oldrtrans0, "oldrtrans");
        _delta_rtrans += clad::getErrorVal(_d_rtrans, _EERepl_rtrans0, "rtrans");
        {
            double _r_d0 = *_d_normr;
            *_d_normr -= _r_d0;
            *_d_normr;
        }
        clad::array<double> _delta_b(_d_b.size());
        int i = 0;
        for (; i < _d_b.size(); i++)
        {
            double _t270 = clad::getErrorVal(_d_b[i], b[i], "b");
            _delta_b[i] += _t270;
            _final_error += _t270;
        }
        i = 0;
        for (; i < _d_x.size(); i++)
        {
            double _t271 = clad::getErrorVal(_d_x[i], _EERepl_x0[i], "x");
            _delta_x[i] += _t271;
            _final_error += _t271;
        }
        double _delta_tolerance = 0;
        _delta_tolerance += clad::getErrorVal(*_d_tolerance, tolerance, "tolerance");
        i = 0;
        for (; i < _d_r.size(); i++)
        {
            double _t272 = clad::getErrorVal(_d_r[i], _EERepl_r0[i], "r");
            _delta_r[i] += _t272;
            _final_error += _t272;
        }
        i = 0;
        for (; i < _d_p.size(); i++)
        {
            double _t273 = clad::getErrorVal(_d_p[i], _EERepl_p0[i], "p");
            _delta_p[i] += _t273;
            _final_error += _t273;
        }
        i = 0;
        for (; i < _d_Ap.size(); i++)
        {
            double _t274 = clad::getErrorVal(_d_Ap[i], _EERepl_Ap0[i], "Ap");
            _delta_Ap[i] += _t274;
            _final_error += _t274;
        }
        clad::array<double> _delta_xexact(_d_xexact.size());
        i = 0;
        for (; i < _d_xexact.size(); i++)
        {
            double _t275 = clad::getErrorVal(_d_xexact[i], xexact[i], "xexact");
            _delta_xexact[i] += _t275;
            _final_error += _t275;
        }
        _final_error += _delta_tolerance + _delta_rtrans + _delta_oldrtrans + _delta_bta + _delta_diff + _delta_alp + _delta_beta + _delta_alpha + _delta_residual;
    }

}