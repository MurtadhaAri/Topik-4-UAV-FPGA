/*
 * run_analysis.c -- System 4 numerical analysis driver.
 *
 *   usage: run_analysis [out_dir=results] [duration_s=60]
 *
 * Produces (every number below is MEASURED by this program, none is typed in):
 *   csv/qformat_table.csv         Q-format properties (closed-form)
 *   csv/cordic_error.csv          CORDIC error vs iteration count
 *   csv/divider_compare.csv       baseline / restoring / non-restoring
 *   csv/quant_hist.csv            rounding/truncation error histogram
 *   csv/pos_integration.csv       sub-LSB position integration experiment
 *   csv/<scenario>_<cfg>_timeseries.csv
 *   csv/summary_rmse.csv  csv/summary_counters.csv  csv/static_criteria.csv
 *   reports/analysis_report.md
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ekf_fixed_adapter.h"
#include "ekf_fixed_divider.h"
#include "ekf_fixed_math.h"
#include "ekf_fixed_trig.h"
#include "s4_compare.h"
#include "s4_dense_ref.h"

static char g_dir[512] = "results";

static FILE *open_out(const char *sub, const char *name)
{
    char path[1024];
    FILE *f;
    snprintf(path, sizeof path, "%s/%s/%s", g_dir, sub, name);
    f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "cannot write %s\n", path);
        exit(2);
    }
    return f;
}

/* ------------------------------------------------------------ Q-format table */
typedef struct { const char *name; int wl, m, n; } qf_t;
static const qf_t k_qf[] = {
    {"Q1.15", 16, 1, 15}, {"Q2.14", 16, 2, 14}, {"Q5.11", 16, 5, 11}, {"Q7.9", 16, 7, 9},
    {"Q1.23", 24, 1, 23}, {"Q2.22", 24, 2, 22}, {"Q7.17", 24, 7, 17}, {"Q8.16", 24, 8, 16},
    {"Q2.30", 32, 2, 30}, {"Q1.31", 32, 1, 31}, {"Q7.25", 32, 7, 25}, {"Q3.29", 32, 3, 29},
    {"Q8.24", 32, 8, 24}, {"Q16.16", 32, 16, 16}};

static void write_qformat_table(void)
{
    FILE *f = open_out("csv", "qformat_table.csv");
    size_t i;
    fprintf(f, "format,word_length,sign_bits,integer_bits_excl_sign,fractional_bits,resolution_lsb,"
               "min_value,max_value,max_rounding_error,max_truncation_error\n");
    for (i = 0; i < sizeof k_qf / sizeof k_qf[0]; i++) {
        const qf_t *q = &k_qf[i];
        double lsb = ldexp(1.0, -q->n);
        fprintf(f, "%s,%d,1,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n", q->name, q->wl, q->m - 1, q->n,
                lsb, -ldexp(1.0, q->m - 1), ldexp(1.0, q->m - 1) - lsb, 0.5 * lsb, lsb);
    }
    fclose(f);
}

/* -------------------------------------------------------------- CORDIC error */
static void write_cordic(void)
{
    static const int its[] = {8, 12, 16, 20, 24, 28, 30};
    FILE *f = open_out("csv", "cordic_error.csv");
    size_t t;
    fprintf(f, "iterations,samples,rms_error_sin,max_error_sin,rms_error_cos,max_error_cos,"
               "max_error_in_lsb_q30\n");
    for (t = 0; t < sizeof its / sizeof its[0]; t++) {
        const int N = 20001;
        double ss = 0, sc = 0, ms = 0, mc = 0;
        int k;
        for (k = 0; k < N; k++) {
            double a = -3.14159265358979 + 2.0 * 3.14159265358979 * k / (N - 1);
            int32_t s, c;
            double es, ec;
            cordic_sincos_n(cordic_angle_from_double(a), (uint8_t)its[t], &s, &c);
            es = fabs(fxp_to_double_q(s, 30) - sin(a));
            ec = fabs(fxp_to_double_q(c, 30) - cos(a));
            ss += es * es;
            sc += ec * ec;
            if (es > ms) ms = es;
            if (ec > mc) mc = ec;
        }
        fprintf(f, "%d,%d,%.6e,%.6e,%.6e,%.6e,%.2f\n", its[t], N, sqrt(ss / N), ms, sqrt(sc / N), mc,
                (ms > mc ? ms : mc) * ldexp(1.0, 30));
    }
    fclose(f);
}

/* -------------------------------------------------------------- divider study */
static uint32_t lcg(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s;
}

static void write_divider(void)
{
    FILE *f = open_out("csv", "divider_compare.csv");
    const div_mode_t modes[3] = {DIV_BASELINE, DIV_RESTORING, DIV_NONRESTORING};
    int m;
    fprintf(f, "mode,samples,mismatch_vs_baseline,max_err_lsb_vs_exact,mean_iterations,"
               "max_iterations,mean_cycles_est,div0_returns_saturated,saturated_results\n");
    for (m = 0; m < 3; m++) {
        uint32_t seed = 12345u;
        const int N = 100000;
        int i, mism = 0, sat0 = 0, satn = 0;
        double maxe = 0, it = 0, cy = 0;
        uint32_t maxit = 0;
        for (i = 0; i < N; i++) {
            int32_t a = (int32_t)(lcg(&seed) >> 8) - (1 << 23); /* +-2^23 raw */
            int32_t b = (int32_t)(lcg(&seed) >> 10) - (1 << 21);
            div_info_t inf;
            int32_t r, rb;
            if (b == 0) b = 1;
            r = div_q(a, b, 16, 16, 16, modes[m], &inf);
            rb = fxp_div_q(a, b, 16, 16, 16);
            if (r != rb) mism++;
            {
                double ex = ((double)a / (double)b) * 65536.0;
                double e = fabs((double)r - ex);
                if (r == INT32_MAX || r == INT32_MIN) satn++;
                else if (e > maxe) maxe = e;
            }
            it += inf.iterations;
            cy += inf.cycles_est;
            if (inf.iterations > maxit) maxit = inf.iterations;
        }
        {
            div_info_t inf;
            int32_t r = div_q(12345, 0, 16, 16, 16, modes[m], &inf);
            sat0 = (inf.div0 && r == INT32_MAX);
        }
        fprintf(f, "%s,%d,%d,%.4f,%.2f,%u,%.2f,%d,%d\n", div_mode_name(modes[m]), N, mism, maxe,
                it / N, maxit, cy / N, sat0, satn);
    }
    fclose(f);
}

/* ------------------------------------------------------- quantization histogram */
static void write_quant_hist(void)
{
    static const struct { const char *name; int frac; } fm[3] = {{"Q1.15", 15}, {"Q1.23", 23}, {"Q2.30", 30}};
    enum { NB = 20 };
    FILE *f = open_out("csv", "quant_hist.csv");
    uint32_t seed = 987654u;
    int fi, i, b;
    fprintf(f, "format,mode,bin_lo_lsb,bin_hi_lsb,count,samples\n");
    for (fi = 0; fi < 3; fi++) {
        long hr[NB], ht[NB];
        const int N = 200000;
        memset(hr, 0, sizeof hr);
        memset(ht, 0, sizeof ht);
        for (i = 0; i < N; i++) {
            double x = ((double)(lcg(&seed) >> 8) / 16777216.0) * 1.9 - 0.95; /* (-0.95,0.95) */
            double lsb = ldexp(1.0, -fm[fi].frac);
            int32_t r = fxp_from_double_q(x, 32u, (uint8_t)fm[fi].frac);
            int32_t t = (int32_t)floor(ldexp(x, fm[fi].frac)); /* truncation (floor) */
            double er = (fxp_to_double_q(r, (uint8_t)fm[fi].frac) - x) / lsb; /* in [-0.5,0.5] */
            double et = (fxp_to_double_q(t, (uint8_t)fm[fi].frac) - x) / lsb; /* in (-1,0]   */
            b = (int)floor((er + 0.5) * NB);
            if (b < 0) b = 0;
            if (b >= NB) b = NB - 1;
            hr[b]++;
            b = (int)floor((et + 1.0) * NB);
            if (b < 0) b = 0;
            if (b >= NB) b = NB - 1;
            ht[b]++;
        }
        for (b = 0; b < NB; b++)
            fprintf(f, "%s,round,%.3f,%.3f,%ld,%d\n", fm[fi].name, -0.5 + (double)b / NB,
                    -0.5 + (double)(b + 1) / NB, hr[b], N);
        for (b = 0; b < NB; b++)
            fprintf(f, "%s,truncate,%.3f,%.3f,%ld,%d\n", fm[fi].name, -1.0 + (double)b / NB,
                    -1.0 + (double)(b + 1) / NB, ht[b], N);
    }
    fclose(f);
}

/* ------------------------- nominal position integration (sub-LSB increments) */
static void write_pos_integration(void)
{
    static const struct { const char *name; int frac; } pf[3] = {{"Q16.16", 16}, {"Q12.20", 20}, {"Q8.24", 24}};
    static const double vs[3] = {0.01, 0.1, 1.0};
    FILE *f = open_out("csv", "pos_integration.csv");
    int fi, vi, k;
    fprintf(f, "pos_format,velocity_mps,steps,dt_s,exact_position_m,fixed_position_m,"
               "abs_error_m,rel_error\n");
    for (fi = 0; fi < 3; fi++)
        for (vi = 0; vi < 3; vi++) {
            const int steps = 30000; /* 60 s @ 500 Hz */
            int32_t dt = fxp_from_double_q(0.002, 32u, 30u);
            int32_t v = fxp_from_double_q(vs[vi], 32u, 24u);
            int32_t p = 0;
            double ex = vs[vi] * ((double)dt / ldexp(1.0, 30)) * steps; /* exact with the same dt */
            for (k = 0; k < steps; k++) {
                int64_t inc = fxp_rshift_round((int64_t)v * dt, 54 - pf[fi].frac);
                p = fxp_sat_wl((int64_t)p + inc, 32u);
            }
            {
                double px = fxp_to_double_q(p, (uint8_t)pf[fi].frac);
                fprintf(f, "%s,%.3f,%d,0.002,%.9f,%.9f,%.9f,%.6e\n", pf[fi].name, vs[vi], steps, ex,
                        px, fabs(px - ex), fabs(px - ex) / ex);
            }
        }
    fclose(f);
}


/* ---------------------------------------------- per-group word-length ablation */
typedef struct { const char *name; const ekf_qcfg_t *base; int dx, U, D, Phi; } variant_t; /* 0 = keep base, else wl */

static void apply_wl(fxp_fmt_t *f, int wl, const ekf_qcfg_t *w16, const ekf_qcfg_t *w24, const ekf_qcfg_t *w32,
                     int which)
{
    const ekf_qcfg_t *src = (wl == 16) ? w16 : (wl == 24) ? w24 : w32;
    switch (which) {
    case 0: *f = src->dx; break;
    case 1: *f = src->U; break;
    case 2: *f = src->D; break;
    default: *f = src->Phi; break;
    }
}

static void write_ablation(double dur, FILE *rp)
{
    const ekf_qcfg_t *w16 = s4_cfg_get(0), *w24 = s4_cfg_get(1), *w32 = s4_cfg_get(2);
    static const variant_t V[] = {
        {"W16 base", NULL, 16, 16, 16, 16},        {"W24 base", NULL, 24, 24, 24, 24},
        {"W32 base", NULL, 32, 32, 32, 32},        {"W24 + D32", NULL, 24, 24, 32, 24},
        {"W24 + D32 + U32", NULL, 24, 32, 32, 24}, {"W32 with D24", NULL, 32, 32, 24, 32},
        {"W32 with U24", NULL, 32, 24, 32, 32},    {"W32 with dx24", NULL, 24, 32, 32, 32},
        {"W32 with Phi24", NULL, 32, 32, 32, 24},  {"W16 + D32", NULL, 16, 16, 32, 16},
        {"mixed-A dx16 U24 D32 Phi16", NULL, 16, 24, 32, 16}, {"mixed-B dx24 U24 D32 Phi24", NULL, 24, 24, 32, 24}};
    static const s4_scenario_t SC[3] = {S4_SC_CONST_VEL, S4_SC_DYNAMIC, S4_SC_HIGH_RANGE};
    FILE *f = open_out("csv", "ablation.csv");
    s4_params_t prm;
    size_t v;
    int s;
    s4_params_default(&prm);
    fprintf(f, "scenario,variant,wl_dx,wl_U,wl_D,wl_Phi,rmse_position_m,rmse_velocity_mps,rmse_attitude_rad,"
               "rmse_P_norm_diag,saturations,D_clamps\n");
    fprintf(rp, "\n## Per-group word-length ablation (fixed vs float RMSE)\n\n"
                "Each row mixes the storage word length of dx / U / D / Phi (16, 24 or 32 bit, formats of the\n"
                "corresponding W16/W24/W32 entry). Scenarios: constvel | dynamic | highrange; each cell: position [m] / attitude [rad].\n\n"
                "| variant | dx | U | D | Phi | constvel pos / att | dynamic pos / att | highrange pos / att |\n"
                "|---|---|---|---|---|---|---|---|\n");
    for (v = 0; v < sizeof V / sizeof V[0]; v++) {
        ekf_qcfg_t c = *w32;
        fprintf(rp, "| %s | %d | %d | %d | %d |", V[v].name, V[v].dx, V[v].U, V[v].D, V[v].Phi);
        c.name = V[v].name;
        apply_wl(&c.dx, V[v].dx, w16, w24, w32, 0);
        apply_wl(&c.U, V[v].U, w16, w24, w32, 1);
        apply_wl(&c.D, V[v].D, w16, w24, w32, 2);
        apply_wl(&c.Phi, V[v].Phi, w16, w24, w32, 3);
        for (s = 0; s < 3; s++) {
            static s4_cmp_t r;
            fxp_counters_t t = {0, 0, 0, 0};
            int g;
            s4_cmp_run_cfg(SC[s], dur, &c, &prm, DIV_BASELINE, NULL, &r);
            for (g = 0; g < (int)FXP_GRP_COUNT; g++) t.sat += r.cnt[g].sat;
            fprintf(f, "%s,%s,%d,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%u,%u\n", s4_scenario_name(SC[s]), V[v].name,
                    V[v].dx, V[v].U, V[v].D, V[v].Phi, s4_acc_rmse(&r.err[0]), s4_acc_rmse(&r.err[1]),
                    s4_acc_rmse(&r.err[2]), s4_acc_rmse(&r.err[12]), t.sat, r.pd_clamp);
            fprintf(rp, " %.2e / %.2e |", s4_acc_rmse(&r.err[0]), s4_acc_rmse(&r.err[2]));
        }
        fprintf(rp, "\n");
    }
    fclose(f);
}

/* ------------------------------------------------------------- scenario runs */
static void write_scenarios(double dur)
{
    FILE *fr = open_out("csv", "summary_rmse.csv");
    FILE *fc = open_out("csv", "summary_counters.csv");
    FILE *fs = open_out("csv", "static_criteria.csv");
    FILE *rp = open_out("reports", "analysis_report.md");
    s4_params_t prm;
    int sc, ci, g;
    static s4_cmp_t res[S4_SC_COUNT][S4_NUM_CFG];

    s4_params_default(&prm);
    fprintf(fr, "scenario,config,kind,group,unit,samples,rmse,mae,max_abs,mean_rel\n");
    fprintf(fc, "scenario,config,var_group,saturation,overflow,underflow,div0\n");
    fprintf(fs, "config,max_position_drift_mm,min_D_norm_float,min_D_norm_fixed,pd_clamp,"
                "d_positive,drift_le_1mm\n");

    for (sc = 0; sc < S4_SC_COUNT; sc++) {
        for (ci = 0; ci < S4_NUM_CFG; ci++) {
            const ekf_qcfg_t *cfg = s4_cfg_get(ci);
            char name[128];
            FILE *ts;
            s4_cmp_t *r = &res[sc][ci];
            snprintf(name, sizeof name, "%s_%s_timeseries.csv", s4_scenario_name((s4_scenario_t)sc),
                     cfg->name);
            ts = open_out("csv", name);
            s4_cmp_run((s4_scenario_t)sc, dur, ci, &prm, DIV_BASELINE, ts, r);
            fclose(ts);
            for (g = 0; g < S4_NGRP; g++) {
                fprintf(fr, "%s,%s,fixed_vs_float,%s,%s,%llu,%.6e,%.6e,%.6e,%.6e\n",
                        s4_scenario_name((s4_scenario_t)sc), cfg->name, S4_GRP_NAME[g], S4_GRP_UNIT[g],
                        (unsigned long long)r->err[g].n, s4_acc_rmse(&r->err[g]), s4_acc_mae(&r->err[g]),
                        s4_acc_max(&r->err[g]), s4_acc_rel(&r->err[g]));
                fprintf(fr, "%s,%s,quantization_roundtrip,%s,%s,%llu,%.6e,%.6e,%.6e,%.6e\n",
                        s4_scenario_name((s4_scenario_t)sc), cfg->name, S4_GRP_NAME[g], S4_GRP_UNIT[g],
                        (unsigned long long)r->qerr[g].n, s4_acc_rmse(&r->qerr[g]),
                        s4_acc_mae(&r->qerr[g]), s4_acc_max(&r->qerr[g]), s4_acc_rel(&r->qerr[g]));
            }
            for (g = 0; g < 3; g++) {
                static const char *tn[3] = {"position", "velocity", "attitude_angle"};
                static const char *tu[3] = {"m", "m/s", "rad"};
                fprintf(fr, "%s,%s,float_sensitivity_floor,%s,%s,%llu,%.6e,%.6e,%.6e,%.6e\n",
                        s4_scenario_name((s4_scenario_t)sc), cfg->name, tn[g], tu[g],
                        (unsigned long long)r->sens[g].n, s4_acc_rmse(&r->sens[g]),
                        s4_acc_mae(&r->sens[g]), s4_acc_max(&r->sens[g]), s4_acc_rel(&r->sens[g]));
                fprintf(fr, "%s,%s,float_vs_truth,%s,%s,%llu,%.6e,%.6e,%.6e,%.6e\n",
                        s4_scenario_name((s4_scenario_t)sc), cfg->name, tn[g], tu[g],
                        (unsigned long long)r->tr_f[g].n, s4_acc_rmse(&r->tr_f[g]),
                        s4_acc_mae(&r->tr_f[g]), s4_acc_max(&r->tr_f[g]), s4_acc_rel(&r->tr_f[g]));
                fprintf(fr, "%s,%s,fixed_vs_truth,%s,%s,%llu,%.6e,%.6e,%.6e,%.6e\n",
                        s4_scenario_name((s4_scenario_t)sc), cfg->name, tn[g], tu[g],
                        (unsigned long long)r->tr_x[g].n, s4_acc_rmse(&r->tr_x[g]),
                        s4_acc_mae(&r->tr_x[g]), s4_acc_max(&r->tr_x[g]), s4_acc_rel(&r->tr_x[g]));
            }
            for (g = 0; g < (int)FXP_GRP_COUNT; g++)
                fprintf(fc, "%s,%s,%s,%u,%u,%u,%u\n", s4_scenario_name((s4_scenario_t)sc), cfg->name,
                        fxp_group_name((fxp_group_t)g), r->cnt[g].sat, r->cnt[g].ovf,
                        r->cnt[g].underflow, r->cnt[g].div0);
            if (sc == S4_SC_STATIC)
                fprintf(fs, "%s,%.6f,%.6e,%.6e,%u,%d,%d\n", cfg->name, r->max_drift_x * 1000.0,
                        r->min_D_f, r->min_D_x, r->pd_clamp, (r->min_D_x > 0.0 && r->pd_clamp == 0),
                        (r->max_drift_x * 1000.0 <= 1.0));
        }
    }
    fclose(fr);
    fclose(fc);
    fclose(fs);
    {
        FILE *fg = open_out("csv", "range_analysis.csv");
        fprintf(fg, "scenario,max_abs_p_m,max_abs_v_mps,max_abs_U_norm_offdiag,max_D_norm,peak_abs_dx_norm,"
                    "int_bits_incl_sign_U,int_bits_incl_sign_D,int_bits_incl_sign_dx\n");
        for (sc = 0; sc < S4_SC_COUNT; sc++) {
            const s4_cmp_t *r0 = &res[sc][2];
            double mu = r0->max_Un, md = r0->max_Dn, mx = r0->peak_dx_n;
            fprintf(fg, "%s,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%d,%d\n", s4_scenario_name((s4_scenario_t)sc),
                    r0->max_abs_p, r0->max_abs_v, mu, md, mx,
                    (mu < 1.0) ? 1 : (int)floor(log2(mu)) + 2, (md < 1.0) ? 1 : (int)floor(log2(md)) + 2,
                    (mx < 1.0) ? 1 : (int)floor(log2(mx)) + 2);
        }
        fclose(fg);
    }

    /* ----------------------------------------------------- markdown report */
    fprintf(rp, "# System 4 -- Analysis Report (auto-generated by tools/run_analysis)\n\n");
    fprintf(rp, "Every number in this file was produced by running the code in this repository on\n"
                "**synthetic** datasets (duration %.0f s, 500 Hz IMU). Nothing is typed in by hand.\n"
                "Model A = float32 study pipeline, Model B = fixed-point study pipeline\n"
                "(PENDING INTEGRATION with the real ES-EKF modules of members 1-3).\n\n", dur);
    fprintf(rp, "Float32 reference accuracy (UDU kernels vs an independent dense double-precision "
                "covariance EKF, worst relative error of diag(P) at the end of the run):\n\n"
                "| scenario | worst relative P_ii error of Model A |\n|---|---|\n");
    for (sc = 0; sc < S4_SC_COUNT; sc++)
        fprintf(rp, "| %s | %.3e |\n", s4_scenario_name((s4_scenario_t)sc),
                s4_verify_float_vs_dense((s4_scenario_t)sc, dur));
    fprintf(rp, "\nFixed-vs-float covariance differences at or below this level cannot be "
                "attributed to the fixed-point model (the reference itself is float32).\n\n");
    fprintf(rp, "Configurations (storage word length of dx / U / D / Phi; nominal path is identical):\n\n"
                "| cfg | dx | U | D | Phi | nominal p / v / q / bias |\n|---|---|---|---|---|---|\n");
    for (ci = 0; ci < S4_NUM_CFG; ci++) {
        const ekf_qcfg_t *c = s4_cfg_get(ci);
        fprintf(rp, "| %s | Q%d.%d | Q%d.%d | Q%d.%d | Q%d.%d | Q%d.%d / Q%d.%d / Q%d.%d / Q%d.%d |\n",
                c->name, c->dx.wl - c->dx.frac, c->dx.frac, c->U.wl - c->U.frac, c->U.frac,
                c->D.wl - c->D.frac, c->D.frac, c->Phi.wl - c->Phi.frac, c->Phi.frac,
                c->p.wl - c->p.frac, c->p.frac, c->v.wl - c->v.frac, c->v.frac, c->q.wl - c->q.frac,
                c->q.frac, c->b.wl - c->b.frac, c->b.frac);
    }
    for (sc = 0; sc < S4_SC_COUNT; sc++) {
        fprintf(rp, "\n## Scenario `%s`\n\nRMSE of fixed vs float reference (log every 0.1 s):\n\n",
                s4_scenario_name((s4_scenario_t)sc));
        fprintf(rp, "| group | unit | W16 RMSE | W24 RMSE | W32 RMSE | W16 max | W24 max | W32 max |\n"
                    "|---|---|---|---|---|---|---|---|\n");
        for (g = 0; g < S4_NGRP; g++)
            fprintf(rp, "| %s | %s | %.3e | %.3e | %.3e | %.3e | %.3e | %.3e |\n", S4_GRP_NAME[g],
                    S4_GRP_UNIT[g], s4_acc_rmse(&res[sc][0].err[g]), s4_acc_rmse(&res[sc][1].err[g]),
                    s4_acc_rmse(&res[sc][2].err[g]), s4_acc_max(&res[sc][0].err[g]),
                    s4_acc_max(&res[sc][1].err[g]), s4_acc_max(&res[sc][2].err[g]));
        fprintf(rp, "\nEstimation error versus TRUTH (RMSE; same sensor data):\n\n"
                    "| quantity | float | W16 | W24 | W32 |\n|---|---|---|---|---|\n");
        {
            static const char *tn[3] = {"position [m]", "velocity [m/s]", "attitude [rad]"};
            for (g = 0; g < 3; g++)
                fprintf(rp, "| %s | %.4e | %.4e | %.4e | %.4e |\n", tn[g],
                        s4_acc_rmse(&res[sc][2].tr_f[g]), s4_acc_rmse(&res[sc][0].tr_x[g]),
                        s4_acc_rmse(&res[sc][1].tr_x[g]), s4_acc_rmse(&res[sc][2].tr_x[g]));
            fprintf(rp, "\nReference sensitivity floor (float vs float with IMU inputs perturbed by "
                        "~1 float32 epsilon, RMSE): position %.3e m, velocity %.3e m/s, attitude %.3e rad. "
                        "A fixed-vs-float difference at or below this level cannot be attributed to "
                        "quantization.\n",
                    s4_acc_rmse(&res[sc][2].sens[0]), s4_acc_rmse(&res[sc][2].sens[1]),
                    s4_acc_rmse(&res[sc][2].sens[2]));
        }
        fprintf(rp, "\nCounters (whole run):\n\n| cfg | saturation | overflow | underflow | div0 | D<=0 clamps | "
                    "Bierman aborts | float D<=0 |\n|---|---|---|---|---|---|---|---|\n");
        for (ci = 0; ci < S4_NUM_CFG; ci++) {
            fxp_counters_t t = {0, 0, 0, 0};
            const s4_cmp_t *r = &res[sc][ci];
            for (g = 0; g < (int)FXP_GRP_COUNT; g++) {
                t.sat += r->cnt[g].sat;
                t.ovf += r->cnt[g].ovf;
                t.underflow += r->cnt[g].underflow;
                t.div0 += r->cnt[g].div0;
            }
            fprintf(rp, "| %s | %u | %u | %u | %u | %u | %u | %u |\n", s4_cfg_get(ci)->name, t.sat,
                    t.ovf, t.underflow, t.div0, r->pd_clamp, r->bierman_abort, r->float_pd_viol);
        }
    }
    fprintf(rp, "\n## ICD static sanity check (Role 5): drift <= 1 mm and D_i > 0\n\n"
                "Applies to the ideal (noise-free) static scenario only. This is a check of the\n"
                "System 4 study pipeline, NOT of the team's full ES-EKF.\n\n"
                "| cfg | max drift [mm] | min D_norm (float) | min D_norm (fixed) | D<=0 clamps | drift<=1mm |\n"
                "|---|---|---|---|---|---|\n");
    for (ci = 0; ci < S4_NUM_CFG; ci++) {
        const s4_cmp_t *r = &res[S4_SC_STATIC][ci];
        fprintf(rp, "| %s | %.6f | %.3e | %.3e | %u | %s |\n", s4_cfg_get(ci)->name,
                r->max_drift_x * 1000.0, r->min_D_f, r->min_D_x, r->pd_clamp,
                (r->max_drift_x * 1000.0 <= 1.0) ? "yes" : "NO");
    }
    {
        FILE *fq = open_out("csv", "quat_norm.csv");
        fprintf(fq, "scenario,config,max_abs_quat_norm_error\n");
        fprintf(rp, "\n## Quaternion normalization (ICD: order-1 Taylor, no sqrt)\n\n"
                    "Maximum | |q| - 1 | of the fixed model over all frames (float reference uses exact sqrt):\n\n"
                    "| scenario | W16 | W24 | W32 |\n|---|---|---|---|\n");
        for (sc = 0; sc < S4_SC_COUNT; sc++) {
            fprintf(rp, "| %s |", s4_scenario_name((s4_scenario_t)sc));
            for (ci = 0; ci < S4_NUM_CFG; ci++) {
                fprintf(fq, "%s,%s,%.6e\n", s4_scenario_name((s4_scenario_t)sc), s4_cfg_get(ci)->name,
                        res[sc][ci].qnorm_dev_x);
                fprintf(rp, " %.3e |", res[sc][ci].qnorm_dev_x);
            }
            fprintf(rp, "\n");
        }
        fclose(fq);
    }
    write_ablation(dur, rp);
    fclose(rp);
}

int main(int argc, char **argv)
{
    double dur = 60.0;
    if (argc > 1) snprintf(g_dir, sizeof g_dir, "%s", argv[1]);
    if (argc > 2) dur = atof(argv[2]);
    if (dur < 1.0 || dur > 60.0) {
        fprintf(stderr, "duration must be within [1,60] s\n");
        return 2;
    }
    printf("[1/6] Q-format table\n");
    write_qformat_table();
    printf("[2/6] CORDIC error\n");
    write_cordic();
    printf("[3/6] divider comparison\n");
    write_divider();
    printf("[4/6] quantization histogram\n");
    write_quant_hist();
    printf("[5/6] position integration experiment\n");
    write_pos_integration();
    printf("[6/6] scenarios (%.0f s each, 4 scenarios x 3 configs)\n", dur);
    write_scenarios(dur);
    printf("done -> %s/{csv,reports}\n", g_dir);
    return 0;
}
