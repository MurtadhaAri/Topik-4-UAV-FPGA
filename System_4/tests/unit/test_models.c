/*
 * Model-level tests of the System 4 STUDY pipelines (Model A float, Model B fixed).
 * Scope: verifies the study kernels and the fixed-vs-float comparison machinery.
 * It does NOT verify the team's ES-EKF (PENDING INTEGRATION).
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "ekf_fixed_adapter.h"
#include "ekf_fixed_divider.h"
#include "s4_compare.h"
#include "s4_dense_ref.h"
#include "test_util.h"

int main(void)
{
    static s4_cmp_t r[S4_SC_COUNT][S4_NUM_CFG];
    static s4_cmp_t again, rr, rn;
    s4_params_t prm;
    int sc, c, g;
    const double DUR = 20.0;
    TEST_BEGIN("study models (float reference / fixed)");

    s4_params_default(&prm);
    CHECK_NEAR(prm.D0[0], 1.0, 0, "ICD D0 position");
    CHECK_NEAR(prm.D0[4], 1e-4, 1e-9, "ICD D0 gyro bias");
    CHECK_NEAR(prm.r_pos, 0.25, 0, "ICD R gnss pos");
    CHECK_NEAR(prm.r_vel, 0.04, 1e-7, "ICD R gnss vel");
    CHECK_NEAR(prm.r_baro, 0.16, 1e-7, "ICD R baro");

    /* 1. verification of the float Thornton/Bierman kernels vs dense double EKF */
    {
        double e1 = s4_verify_float_vs_dense(S4_SC_STATIC, 30.0);
        double e2 = s4_verify_float_vs_dense(S4_SC_DYNAMIC, 30.0);
        printf("  float UDU vs dense double P: worst relative diag error static %.2e dynamic %.2e\n", e1, e2);
        CHECK(e1 < 1e-3, "static: UDU kernels match dense covariance EKF");
        CHECK(e2 < 1e-3, "dynamic: UDU kernels match dense covariance EKF");
    }

    /* 2. run everything */
    for (sc = 0; sc < S4_SC_COUNT; sc++)
        for (c = 0; c < S4_NUM_CFG; c++) {
            CHECK_EQ(s4_cmp_run((s4_scenario_t)sc, DUR, c, &prm, DIV_BASELINE, NULL, &r[sc][c]), 0, "run ok");
            CHECK(!r[sc][c].nan_flag, "no NaN/Inf in states");
            CHECK_EQ(r[sc][c].float_pd_viol, 0, "float reference keeps D>0");
            CHECK_EQ(r[sc][c].bierman_abort, 0, "no Bierman aborts");
            CHECK(r[sc][c].n_log == (int)(DUR / 0.1), "log tick count");
        }

    /* 3. ICD static sanity (study pipeline only): drift <= 1 mm, D_i > 0 */
    for (c = 0; c < S4_NUM_CFG; c++) {
        CHECK(r[S4_SC_STATIC][c].max_drift_f * 1000.0 <= 1.0, "static: float drift <= 1 mm");
        CHECK(r[S4_SC_STATIC][c].max_drift_x * 1000.0 <= 1.0, "static: fixed drift <= 1 mm");
        CHECK(r[S4_SC_STATIC][c].min_D_x > 0.0, "static: fixed D_i > 0");
        CHECK_EQ(r[S4_SC_STATIC][c].pd_clamp, 0, "static: no positive-definiteness clamp");
    }

    /* 4. word-length ordering of covariance fidelity (measured property, 3 scenarios) */
    for (sc = 0; sc < S4_SC_COUNT; sc++) {
        double e16 = s4_acc_rmse(&r[sc][0].err[12]);
        double e24 = s4_acc_rmse(&r[sc][1].err[12]);
        double e32 = s4_acc_rmse(&r[sc][2].err[12]);
        printf("  %-9s P_norm_diag RMSE  W16 %.3e  W24 %.3e  W32 %.3e\n", s4_scenario_name((s4_scenario_t)sc), e16, e24, e32);
        CHECK(e16 > e24, "P diag error: W16 > W24");
        CHECK(e32 < 1e-3, "W32 P diag RMSE < 1e-3");
    }

    /* 4b. ablation: D storage width is the dominant accuracy lever (measured, robust finding).
     *     W24 with only D widened to 32 bit must be >= 2x better in position than plain W24. */
    for (sc = S4_SC_CONST_VEL; sc <= S4_SC_HIGH_RANGE; sc++) {
        ekf_qcfg_t mix = *s4_cfg_get(1);
        static s4_cmp_t rm;
        mix.D = s4_cfg_get(2)->D;
        mix.name = "W24+D32";
        s4_cmp_run_cfg((s4_scenario_t)sc, DUR, &mix, &prm, DIV_BASELINE, NULL, &rm);
        printf("  %-9s position RMSE  W24 %.2e  W24+D32 %.2e  (%.1fx better)\n", s4_scenario_name((s4_scenario_t)sc),
               s4_acc_rmse(&r[sc][1].err[0]), s4_acc_rmse(&rm.err[0]),
               s4_acc_rmse(&r[sc][1].err[0]) / s4_acc_rmse(&rm.err[0]));
        CHECK(s4_acc_rmse(&rm.err[0]) * 2.0 < s4_acc_rmse(&r[sc][1].err[0]), "widening only D to 32 bit halves W24 position error");
    }

    /* 5. W32 acceptance on the realistic scenarios (PROPOSED thresholds) */
    for (sc = S4_SC_CONST_VEL; sc <= S4_SC_DYNAMIC; sc++) {
        const s4_cmp_t *x = &r[sc][2];
        printf("  %-9s W32 vs float: pos %.2e m  vel %.2e m/s  att %.2e rad\n", s4_scenario_name((s4_scenario_t)sc),
               s4_acc_rmse(&x->err[0]), s4_acc_rmse(&x->err[1]), s4_acc_rmse(&x->err[2]));
        CHECK(s4_acc_rmse(&x->err[0]) < 1e-2, "W32 position RMSE vs float < 1 cm");
        CHECK(s4_acc_rmse(&x->err[1]) < 1e-2, "W32 velocity RMSE vs float < 1 cm/s");
        CHECK(s4_acc_rmse(&x->err[2]) < 1e-3, "W32 attitude RMSE vs float < 1 mrad");
    }
    /* highrange: the float reference is quiet (sensitivity floor) and W32 saturates nothing */
    CHECK(s4_acc_rmse(&r[S4_SC_HIGH_RANGE][2].sens[0]) < 1e-3, "float reference sensitivity floor is tiny");
    {
        fxp_counters_t t = {0, 0, 0, 0};
        for (g = 0; g < (int)FXP_GRP_COUNT; g++) t.sat += r[S4_SC_HIGH_RANGE][2].cnt[g].sat;
        CHECK_EQ(t.sat, 0, "highrange W32: no saturation with Q4.x dx headroom");
    }
    /* Newton inverse-sqrt normalization keeps |q| = 1 even after >0.5 rad injections
     * (the ICD order-1 Taylor rule, -DS4_QNORM_TAYLOR, gave 0.44 on highrange). */
    for (sc = S4_SC_CONST_VEL; sc <= S4_SC_HIGH_RANGE; sc++)
        CHECK(r[sc][2].qnorm_dev_x < 1e-5, "quaternion norm within 1e-5 in all scenarios");
    printf("  highrange max ||q|-1| (W32) = %.3e\n", r[S4_SC_HIGH_RANGE][2].qnorm_dev_x);
    /* the float estimate is a sane estimator (so comparisons mean something) */
    CHECK(s4_acc_rmse(&r[S4_SC_CONST_VEL][2].tr_f[0]) < 0.5, "constvel float position RMSE vs truth < 0.5 m");
    CHECK(s4_acc_rmse(&r[S4_SC_DYNAMIC][2].tr_f[0]) < 0.5, "dynamic float position RMSE vs truth < 0.5 m");

    /* 6. determinism and divider-structure invariance (bit exact) */
    s4_cmp_run(S4_SC_DYNAMIC, DUR, 1, &prm, DIV_BASELINE, NULL, &again);
    CHECK(s4_acc_rmse(&again.err[0]) == s4_acc_rmse(&r[S4_SC_DYNAMIC][1].err[0]), "deterministic re-run");
    s4_cmp_run(S4_SC_DYNAMIC, DUR, 1, &prm, DIV_RESTORING, NULL, &rr);
    s4_cmp_run(S4_SC_DYNAMIC, DUR, 1, &prm, DIV_NONRESTORING, NULL, &rn);
    for (g = 0; g < S4_NGRP; g++) {
        CHECK(s4_acc_rmse(&rr.err[g]) == s4_acc_rmse(&again.err[g]), "restoring divider == baseline results");
        CHECK(s4_acc_rmse(&rn.err[g]) == s4_acc_rmse(&again.err[g]), "non-restoring divider == baseline results");
    }
    CHECK_EQ(div_get_pipeline_mode(), DIV_BASELINE, "pipeline divider mode restored");
    return TEST_END();
}
