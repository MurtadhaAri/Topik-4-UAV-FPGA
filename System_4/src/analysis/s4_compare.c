/*
 * s4_compare.c -- lock-step comparison of the float and fixed study models.
 */
#include "s4_compare.h"

#include <math.h>
#include <string.h>

#include "ekf_fixed_adapter.h"

const char *const S4_GRP_NAME[S4_NGRP] = {"position", "velocity", "attitude_angle", "accel_bias",
                                           "gyro_bias", "quaternion",  "U_norm",         "D_norm_pos",
                                           "D_norm_vel", "D_norm_att", "D_norm_ba",       "D_norm_bg",
                                           "P_norm_diag", "P_norm_offdiag"};
const char *const S4_GRP_UNIT[S4_NGRP] = {"m", "m/s", "rad", "m/s^2", "rad/s", "-", "-",
                                          "-", "-",   "-",   "-",     "-", "-", "-"};

static sensor_frame_t g_frames[S4_MAX_FRAMES];
static s4_truth_t g_truth[S4_MAX_FRAMES];

typedef struct {
    double p[3], v[3], q[4], ba[3], bg[3];
    double Dn[ERR_DIM];
    double Un[U_UPPER_SIZE];
    double Pn[ERR_DIM][ERR_DIM]; /* normalized covariance U D U^T */
} dvec_t;

static void dvec_cov(dvec_t *d)
{
    int i, j, k;
    for (i = 0; i < ERR_DIM; i++)
        for (j = i; j < ERR_DIM; j++) {
            double s = 0.0;
            for (k = j; k < ERR_DIM; k++) s += d->Un[U_IDX(i, k)] * d->Dn[k] * d->Un[U_IDX(j, k)];
            d->Pn[i][j] = d->Pn[j][i] = s;
        }
}

static void dvec_float(const es_ekf_t *e, dvec_t *d)
{
    int i, j;
    for (i = 0; i < 3; i++) {
        d->p[i] = e->x_nom.p[i];
        d->v[i] = e->x_nom.v[i];
        d->ba[i] = e->x_nom.ba[i];
        d->bg[i] = e->x_nom.bg[i];
    }
    for (i = 0; i < 4; i++) d->q[i] = e->x_nom.q[i];
    for (i = 0; i < ERR_DIM; i++) d->Dn[i] = s4_float_Dn(e, i);
    for (j = 0; j < ERR_DIM; j++)
        for (i = 0; i <= j; i++) d->Un[U_IDX(i, j)] = s4_float_Un(e, i, j);
    dvec_cov(d);
}

static void dvec_fixed(const es_ekf_fixed_t *e, const ekf_qcfg_t *c, dvec_t *d)
{
    int i, j;
    for (i = 0; i < 3; i++) {
        d->p[i] = fxp_to_double_q(e->x_nom.p[i], c->p.frac);
        d->v[i] = fxp_to_double_q(e->x_nom.v[i], c->v.frac);
        d->ba[i] = fxp_to_double_q(e->x_nom.ba[i], c->b.frac);
        d->bg[i] = fxp_to_double_q(e->x_nom.bg[i], c->b.frac);
    }
    for (i = 0; i < 4; i++) d->q[i] = fxp_to_double_q(e->x_nom.q[i], c->q.frac);
    for (i = 0; i < ERR_DIM; i++) d->Dn[i] = s4_fixed_Dn(e, c, i);
    for (j = 0; j < ERR_DIM; j++)
        for (i = 0; i <= j; i++) d->Un[U_IDX(i, j)] = s4_fixed_Un(e, c, i, j);
    dvec_cov(d);
}

static void accumulate(s4_acc_t *acc, const dvec_t *ref, const dvec_t *tst)
{
    int i, j, sg;
    for (i = 0; i < 3; i++) {
        s4_acc_add(&acc[0], ref->p[i], tst->p[i]);
        s4_acc_add(&acc[1], ref->v[i], tst->v[i]);
        s4_acc_add(&acc[3], ref->ba[i], tst->ba[i]);
        s4_acc_add(&acc[4], ref->bg[i], tst->bg[i]);
    }
    s4_acc_add(&acc[2], 0.0, s4_quat_angle(ref->q, tst->q)); /* rotation angle difference */
    sg = s4_quat_sign(ref->q, tst->q);
    for (i = 0; i < 4; i++) s4_acc_add(&acc[5], ref->q[i], sg * tst->q[i]);
    for (j = 1; j < ERR_DIM; j++)
        for (i = 0; i < j; i++) s4_acc_add(&acc[6], ref->Un[U_IDX(i, j)], tst->Un[U_IDX(i, j)]);
    for (i = 0; i < ERR_DIM; i++) s4_acc_add(&acc[7 + i / 3], ref->Dn[i], tst->Dn[i]);
    for (i = 0; i < ERR_DIM; i++)
        for (j = i; j < ERR_DIM; j++)
            s4_acc_add(&acc[(i == j) ? 12 : 13], ref->Pn[i][j], tst->Pn[i][j]);
}

static int dvec_finite(const dvec_t *d)
{
    int i;
    for (i = 0; i < 3; i++)
        if (!isfinite(d->p[i]) || !isfinite(d->v[i])) return 0;
    for (i = 0; i < 4; i++)
        if (!isfinite(d->q[i])) return 0;
    return 1;
}

int s4_cmp_run(s4_scenario_t sc, double duration_s, int cfg_idx, const s4_params_t *prm,
               div_mode_t divmode, FILE *ts_csv, s4_cmp_t *out)
{
    const ekf_qcfg_t *cfg = s4_cfg_get(cfg_idx);
    if (!cfg) return -1;
    return s4_cmp_run_cfg(sc, duration_s, cfg, prm, divmode, ts_csv, out);
}

int s4_cmp_run_cfg(s4_scenario_t sc, double duration_s, const ekf_qcfg_t *cfg,
                   const s4_params_t *prm, div_mode_t divmode, FILE *ts_csv, s4_cmp_t *out)
{
    static s4f_model_t mf;
    static s4f_model_t mf2; /* perturbed-input twin of the float reference */
    static s4x_model_t mx;
    uint32_t lcg = 2463534242u;
    int n, k, g, i, j;
    div_mode_t old_mode = div_get_pipeline_mode();

    memset(out, 0, sizeof *out);
    for (g = 0; g < S4_NGRP; g++) {
        s4_acc_init(&out->err[g]);
        s4_acc_init(&out->qerr[g]);
    }
    for (g = 0; g < 3; g++) {
        s4_acc_init(&out->sens[g]);
        s4_acc_init(&out->tr_f[g]);
        s4_acc_init(&out->tr_x[g]);
    }
    out->min_D_f = out->min_D_x = 1e300;

    n = s4_dataset_generate(sc, duration_s, g_frames, g_truth, S4_MAX_FRAMES);
    out->n_frames = n;

    div_set_pipeline_mode(divmode);
    fxp_stats_reset();
    s4f_init(&mf, prm);
    s4x_init(&mx, cfg, prm);
    {
        float v0[3];
        s4_dataset_filter_init(sc, v0);
        s4f_init(&mf2, prm);
        s4f_set_velocity(&mf2, v0);
        s4f_set_velocity(&mf, v0);
        s4x_set_velocity(&mx, v0);
    }

    if (ts_csv) {
        fprintf(ts_csv,
                "t,truth_pN,truth_pE,truth_pD,fl_pN,fl_pE,fl_pD,fx_pN,fx_pE,fx_pD,"
                "truth_vN,truth_vE,truth_vD,fl_vN,fl_vE,fl_vD,fx_vN,fx_vE,fx_vD,"
                "fl_roll,fl_pitch,fl_yaw,fx_roll,fx_pitch,fx_yaw,"
                "err_pos,err_vel,err_att,err_ba,err_bg,sat_total,ovf_total,max_rel_D\n");
    }

    for (k = 0; k < n; k++) {
        sensor_frame_fixed_t ff;
        double dpf, dpx;
        s4f_step(&mf, &g_frames[k]);
        {
            sensor_frame_t pf = g_frames[k];
            for (i = 0; i < 3; i++) {
                lcg = lcg * 1664525u + 1013904223u;
                pf.imu_acc[i] *= 1.0f + 1.1920929e-7f * ((float)(lcg >> 8) / 8388608.0f - 1.0f);
                lcg = lcg * 1664525u + 1013904223u;
                pf.imu_gyr[i] *= 1.0f + 1.1920929e-7f * ((float)(lcg >> 8) / 8388608.0f - 1.0f);
            }
            s4f_step(&mf2, &pf);
        }
        sensor_frame_to_fixed(&g_frames[k], &ff, cfg);
        s4x_step(&mx, &ff);

        {
            double s = 0.0, dev;
            for (i = 0; i < 4; i++) {
                double q = fxp_to_double_q(mx.ekf.x_nom.q[i], cfg->q.frac);
                s += q * q;
            }
            dev = fabs(sqrt(s) - 1.0);
            if (dev > out->qnorm_dev_x) out->qnorm_dev_x = dev;
        }
        dpf = dpx = 0.0;
        for (i = 0; i < 3; i++) {
            double a = fabs((double)mf.ekf.x_nom.p[i]);
            double b = fabs(fxp_to_double_q(mx.ekf.x_nom.p[i], cfg->p.frac));
            if (a > dpf) dpf = a;
            if (b > dpx) dpx = b;
        }
        if (dpf > out->max_drift_f) out->max_drift_f = dpf;
        if (dpx > out->max_drift_x) out->max_drift_x = dpx;

        if ((k + 1) % S4_LOG_EVERY == 0) {
            dvec_t ref, tst, rt;
            es_ekf_fixed_t tmp;
            s4_acc_t tick[S4_NGRP];
            double maxrel = 0.0, rpy_f[3], rpy_x[3];
            fxp_counters_t tot;

            dvec_float(&mf.ekf, &ref);
            dvec_fixed(&mx.ekf, cfg, &tst);
            if (!dvec_finite(&ref) || !dvec_finite(&tst)) out->nan_flag = 1;
            accumulate(out->err, &ref, &tst);
            ekf_float_to_fixed(&mf.ekf, &tmp, cfg);
            dvec_fixed(&tmp, cfg, &rt);
            accumulate(out->qerr, &ref, &rt);
            out->n_log++;
            for (i = 0; i < 3; i++) {
                s4_acc_add(&out->tr_f[0], g_truth[k].p[i], ref.p[i]);
                s4_acc_add(&out->tr_x[0], g_truth[k].p[i], tst.p[i]);
                s4_acc_add(&out->tr_f[1], g_truth[k].v[i], ref.v[i]);
                s4_acc_add(&out->tr_x[1], g_truth[k].v[i], tst.v[i]);
            }
            for (i = 0; i < 3; i++) {
                s4_acc_add(&out->sens[0], ref.p[i], mf2.ekf.x_nom.p[i]);
                s4_acc_add(&out->sens[1], ref.v[i], mf2.ekf.x_nom.v[i]);
            }
            {
                double q2[4];
                for (i = 0; i < 4; i++) q2[i] = mf2.ekf.x_nom.q[i];
                s4_acc_add(&out->sens[2], 0.0, s4_quat_angle(ref.q, q2));
            }
            s4_acc_add(&out->tr_f[2], 0.0, s4_quat_angle(g_truth[k].q, ref.q));
            s4_acc_add(&out->tr_x[2], 0.0, s4_quat_angle(g_truth[k].q, tst.q));

            for (j = 1; j < ERR_DIM; j++)
                for (i = 0; i < j; i++)
                    if (fabs(ref.Un[U_IDX(i, j)]) > out->max_Un) out->max_Un = fabs(ref.Un[U_IDX(i, j)]);
            for (i = 0; i < ERR_DIM; i++)
                if (ref.Dn[i] > out->max_Dn) out->max_Dn = ref.Dn[i];
            for (i = 0; i < 3; i++) {
                if (fabs(ref.p[i]) > out->max_abs_p) out->max_abs_p = fabs(ref.p[i]);
                if (fabs(ref.v[i]) > out->max_abs_v) out->max_abs_v = fabs(ref.v[i]);
            }
            for (i = 0; i < ERR_DIM; i++) {
                double rel = (ref.Dn[i] > 0.0) ? fabs(tst.Dn[i] - ref.Dn[i]) / ref.Dn[i] : 0.0;
                if (rel > maxrel) maxrel = rel;
                if (ref.Dn[i] < out->min_D_f) out->min_D_f = ref.Dn[i];
                if (tst.Dn[i] < out->min_D_x) out->min_D_x = tst.Dn[i];
            }

            if (ts_csv) {
                for (g = 0; g < S4_NGRP; g++) s4_acc_init(&tick[g]);
                accumulate(tick, &ref, &tst);
                tot = fxp_stats_total();
                s4_quat_to_euler(ref.q, rpy_f);
                s4_quat_to_euler(tst.q, rpy_x);
                fprintf(ts_csv, "%.3f", (k + 1) * S4_IMU_DT);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", g_truth[k].p[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", ref.p[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", tst.p[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", g_truth[k].v[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", ref.v[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", tst.v[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", rpy_f[i]);
                for (i = 0; i < 3; i++) fprintf(ts_csv, ",%.9g", rpy_x[i]);
                /* per-tick error magnitudes (Euclidean norm / angle) */
                {
                    double ep = 0, ev = 0, eba = 0, ebg = 0;
                    for (i = 0; i < 3; i++) {
                        ep += (tst.p[i] - ref.p[i]) * (tst.p[i] - ref.p[i]);
                        ev += (tst.v[i] - ref.v[i]) * (tst.v[i] - ref.v[i]);
                        eba += (tst.ba[i] - ref.ba[i]) * (tst.ba[i] - ref.ba[i]);
                        ebg += (tst.bg[i] - ref.bg[i]) * (tst.bg[i] - ref.bg[i]);
                    }
                    fprintf(ts_csv, ",%.9g,%.9g,%.9g,%.9g,%.9g,%u,%u,%.9g\n", sqrt(ep), sqrt(ev),
                            s4_quat_angle(ref.q, tst.q), sqrt(eba), sqrt(ebg), tot.sat, tot.ovf,
                            maxrel);
                }
            }
        }
    }

    for (g = 0; g < (int)FXP_GRP_COUNT; g++) out->cnt[g] = fxp_stats_get((fxp_group_t)g);
    out->peak_dx_n = mf.peak_dx_n;
    out->pd_clamp = mx.n_pd_clamp;
    out->bierman_abort = mx.n_bierman_abort;
    out->float_pd_viol = mf.n_pd_viol;
    div_set_pipeline_mode(old_mode);
    return 0;
}
