/* Adapter + type-contract tests (ICD constants, layout, float<->fixed round trips). */
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "ekf_fixed_adapter.h"
#include "ekf_fixed_math.h"
#include "test_util.h"

static double rnd(double lo, double hi) { return lo + (hi - lo) * ((double)(t_rand() >> 8) / 16777216.0); }

int main(void)
{
    int c, i, j, rep;
    TEST_BEGIN("ekf_fixed_adapter / ICD contract");

    /* ICD contract (Bagian I, butir 5) */
    CHECK_EQ(NOM_DIM, 16, "NOM_DIM");
    CHECK_EQ(ERR_DIM, 15, "ERR_DIM");
    CHECK_EQ(U_UPPER_SIZE, 120, "U_UPPER_SIZE");
    CHECK_EQ(U_UPPER_SIZE, ERR_DIM * (ERR_DIM + 1) / 2, "U_UPPER_SIZE formula");
    CHECK_EQ(U_IDX(0, 0), 0, "U_IDX(0,0)");
    CHECK_EQ(U_IDX(0, 1), 1, "U_IDX(0,1)");
    CHECK_EQ(U_IDX(1, 1), 2, "U_IDX(1,1)");
    CHECK_EQ(U_IDX(14, 14), 119, "U_IDX(14,14) = last element");
    {
        int seen[U_UPPER_SIZE];
        memset(seen, 0, sizeof seen);
        for (j = 0; j < ERR_DIM; j++)
            for (i = 0; i <= j; i++) seen[U_IDX(i, j)]++;
        for (i = 0; i < U_UPPER_SIZE; i++)
            if (seen[i] != 1) { CHECK(0, "U_IDX is not a bijection onto 0..119"); break; }
        CHECK(1, "U_IDX bijection");
    }
    CHECK_EQ(sizeof(((es_ekf_fixed_t *)0)->U) / sizeof(int32_t), U_UPPER_SIZE, "fixed U has 120 elements");
    CHECK_EQ(sizeof(((es_ekf_fixed_t *)0)->D) / sizeof(int32_t), ERR_DIM, "fixed D has 15 elements");
    CHECK_EQ(sizeof(((es_ekf_fixed_t *)0)->dx) / sizeof(int32_t), ERR_DIM, "fixed dx has 15 elements");
    CHECK_EQ(sizeof(nominal_state_fixed_t) / sizeof(int32_t), NOM_DIM, "nominal fixed = 16 words");
    CHECK_EQ(sizeof(nominal_state_t) / sizeof(float), NOM_DIM, "nominal float = 16 words");
    CHECK_EQ(offsetof(nominal_state_fixed_t, q), 6 * sizeof(int32_t), "q starts at index 6");
    CHECK_EQ(offsetof(nominal_state_fixed_t, ba), 10 * sizeof(int32_t), "ba starts at index 10");
    CHECK_EQ(offsetof(nominal_state_fixed_t, bg), 13 * sizeof(int32_t), "bg starts at index 13");
    CHECK_NEAR(S4_GRAVITY_MPS2, 9.80665, 0.0, "ICD gravity constant");

    /* scale factors cover the ICD bounds */
    CHECK(s4_scale(0) >= 10.0 && s4_scale(3) >= 10.0, "pos/vel scale >= 10");
    CHECK(s4_scale(6) >= 0.5, "attitude scale >= 0.5 rad");
    CHECK(s4_scale(9) >= 2.0, "accel-bias scale >= 2");
    CHECK(s4_scale(12) >= 0.2, "gyro-bias scale >= 0.2");
    CHECK_NEAR(s4_scale(0) / 10.0, 1.6, 1e-12, "pow2 scale costs 1.6x headroom on position");

    for (c = 0; c < S4_NUM_CFG; c++) {
        const ekf_qcfg_t *cfg = s4_cfg_get(c);
        double lsb_dx = ldexp(1.0, -cfg->dx.frac), lsb_U = ldexp(1.0, -cfg->U.frac);
        double lsb_D = ldexp(1.0, -cfg->D.frac);
        double mp = 0, mdx = 0, mU = 0, mD = 0, mq = 0, mv = 0;
        printf("  cfg %s\n", cfg->name);
        CHECK(cfg->dx.wl == cfg->U.wl && cfg->U.wl == cfg->D.wl, "same word length on the error path");
        for (rep = 0; rep < 200; rep++) {
            es_ekf_t a, b;
            es_ekf_fixed_t fx;
            double e;
            memset(&a, 0, sizeof a);
            for (i = 0; i < 3; i++) {
                a.x_nom.p[i] = (float)rnd(-1000, 1000);
                a.x_nom.v[i] = (float)rnd(-50, 50);
                a.x_nom.ba[i] = (float)rnd(-2, 2);
                a.x_nom.bg[i] = (float)rnd(-0.2, 0.2);
            }
            {
                double q[4], n = 0;
                for (i = 0; i < 4; i++) { q[i] = rnd(-1, 1); n += q[i] * q[i]; }
                n = sqrt(n);
                for (i = 0; i < 4; i++) a.x_nom.q[i] = (float)(q[i] / n);
            }
            for (i = 0; i < ERR_DIM; i++) {
                a.dx[i] = (float)(rnd(-1, 1) * s4_scale(i));
                a.D[i] = (float)(rnd(0.001, 0.5) * s4_scale(i) * s4_scale(i));
            }
            for (j = 0; j < ERR_DIM; j++)
                for (i = 0; i <= j; i++)
                    a.U[U_IDX(i, j)] = (i == j) ? 1.0f : (float)(rnd(-8, 8) * s4_scale(i) / s4_scale(j));
            ekf_float_to_fixed(&a, &fx, cfg);
            ekf_fixed_to_float(&fx, &b, cfg);
            for (i = 0; i < 3; i++) {
                e = fabs((double)a.x_nom.p[i] - b.x_nom.p[i]); if (e > mp) mp = e;
                e = fabs((double)a.x_nom.v[i] - b.x_nom.v[i]); if (e > mv) mv = e;
            }
            for (i = 0; i < 4; i++) { e = fabs((double)a.x_nom.q[i] - b.x_nom.q[i]); if (e > mq) mq = e; }
            for (i = 0; i < ERR_DIM; i++) {
                e = fabs(ldexp((double)a.dx[i] - b.dx[i], -s4_scale_log2[i])); if (e > mdx) mdx = e;
                e = fabs(ldexp((double)a.D[i] - b.D[i], -2 * s4_scale_log2[i])); if (e > mD) mD = e;
            }
            for (j = 0; j < ERR_DIM; j++)
                for (i = 0; i <= j; i++) {
                    e = fabs(ldexp((double)a.U[U_IDX(i, j)] - b.U[U_IDX(i, j)], s4_scale_log2[j] - s4_scale_log2[i]));
                    if (e > mU) mU = e;
                }
            for (j = 0; j < ERR_DIM; j++)
                if (fx.U[U_IDX(j, j)] != (int32_t)(1 << cfg->U.frac)) { CHECK(0, "U diagonal is not exactly 1.0"); break; }
        }
        printf("    max round-trip error: p %.3e v %.3e q %.3e dx_n %.3e (lsb %.3e) U_n %.3e (lsb %.3e) D_n %.3e (lsb %.3e)\n",
               mp, mv, mq, mdx, lsb_dx, mU, lsb_U, mD, lsb_D);
        /* float32 itself contributes ~6e-8 relative; allow it on top of the half-LSB bound */
        CHECK(mp <= 0.5 * ldexp(1.0, -cfg->p.frac) + 1000 * 6e-8, "position within 1/2 LSB (+float32 eps)");
        CHECK(mq <= 0.5 * ldexp(1.0, -cfg->q.frac) + 1.5e-7, "quaternion within 1/2 LSB (+float32 eps)");
        CHECK(mdx <= 0.5 * lsb_dx + 1.5e-7, "dx_n within 1/2 LSB (+float32 eps)");
        CHECK(mU <= 0.5 * lsb_U + 8 * 6e-8 * 8, "U_n within 1/2 LSB (+float32 eps)");
        CHECK(mD <= 0.5 * lsb_D + 1.5e-7, "D_n within 1/2 LSB (+float32 eps)");
    }

    /* quaternion +1 does not wrap in Q2.30 */
    {
        es_ekf_t a;
        es_ekf_fixed_t fx;
        memset(&a, 0, sizeof a);
        a.x_nom.q[0] = 1.0f;
        for (i = 0; i < ERR_DIM; i++) { a.U[U_IDX(i, i)] = 1.0f; a.D[i] = 0.1f; }
        ekf_float_to_fixed(&a, &fx, s4_cfg_get(2));
        CHECK_EQ(fx.x_nom.q[0], 1 << 30, "qw=+1 exactly representable in Q2.30");
        CHECK(fx.x_nom.q[0] > 0, "no wrap-around");
    }

    /* sensor frame conversion */
    {
        sensor_frame_t f;
        sensor_frame_fixed_t ff;
        const ekf_qcfg_t *cfg = s4_cfg_get(1);
        memset(&f, 0, sizeof f);
        f.dt = 0.002f;
        f.type = SENS_GNSS_POS_VEL;
        f.imu_acc[2] = -9.80665f;
        f.imu_gyr[0] = 0.5f;
        f.meas[0] = 123.456f;
        f.meas[3] = -2.5f;
        sensor_frame_to_fixed(&f, &ff, cfg);
        CHECK_NEAR(fxp_to_double_q(ff.dt, 30), 0.002, 1e-9, "dt Q2.30");
        CHECK_NEAR(fxp_to_double_q(ff.imu_acc[2], 24), -9.80665, 1e-6, "acc Q8.24");
        CHECK_NEAR(fxp_to_double_q(ff.imu_gyr[0], 24), 0.5, 1e-7, "gyr Q8.24");
        CHECK_NEAR(fxp_to_double_q(ff.meas[0], cfg->p.frac), 123.456, 1e-4, "meas pos in p format");
        CHECK_NEAR(fxp_to_double_q(ff.meas[3], cfg->v.frac), -2.5, 1e-7, "meas vel in v format");
        CHECK_EQ(ff.type, SENS_GNSS_POS_VEL, "type preserved");
    }
    return TEST_END();
}
