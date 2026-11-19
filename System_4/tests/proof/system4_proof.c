/*
 * system4_proof.c -- ONE-FILE DUMMY-TEST PROOF of the System 4 best configuration.
 *
 *   usage: system4_proof [results_dir=results]      exit code 0 = all proofs hold
 *
 * What it proves (everything is computed here, nothing is typed in):
 *   A. dummy arithmetic checks   saturation, rounding, +1 in Q2.30, CORDIC, divider
 *   B. best-configuration search every candidate (16/24/32-bit mixes) is run through
 *                                the same DUT harness on 4 synthetic scenarios; the
 *                                cheapest candidate meeting ALL acceptance thresholds
 *                                (system4.h) is selected
 *   C. header consistency        the selection must equal s4_best_cfg() in system4.h
 *   D. negative control          the TODO algorithm template (S4_DUT_TEMPLATE) must be
 *                                REJECTED by the same harness -> the harness can fail
 *   E. report                    <results>/proof/best_config_proof.{md,csv}
 *
 * Scope: study pipeline vs float32 reference on SYNTHETIC data.  It does not prove the
 * team's ES-EKF (PENDING INTEGRATION).  Thresholds are PROPOSED (system4.h).
 *
 * TEMPLATE USE: write your DUT as an s4_dut_t (see S4_DUT_TEMPLATE in
 * src/model/system4_dut.c), then replace S4_DUT_STUDY in run_candidates().
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "s4_dataset.h"
#include "s4_metrics.h"
#include "system4.h"

#define DUR_S 60.0

static int g_checks = 0, g_fails = 0;
#define PROOF(cond, msg)                                                          \
    do {                                                                          \
        g_checks++;                                                               \
        if (!(cond)) {                                                            \
            g_fails++;                                                            \
            printf("  [FAIL] %s\n", (msg));                                       \
        }                                                                         \
    } while (0)

/* ----------------------------------------------------------- DUT harness */
typedef struct {
    double rmse_p, rmse_v, rmse_a, qdev, drift_mm;
    uint32_t sat, ovf, pd;
    int nan_flag, d_positive;
} metrics_t;

static sensor_frame_t g_frames[S4_MAX_FRAMES];
static s4_truth_t g_truth[S4_MAX_FRAMES];
static union {
    unsigned char raw[S4_DUT_CTX_MAX];
    double align_d;
    long long align_l;
} g_ctx;

static metrics_t run_dut(const s4_dut_t *dut, const ekf_qcfg_t *cfg, s4_scenario_t sc)
{
    static s4f_model_t mf;
    s4_params_t prm;
    s4_acc_t ap, av, aa;
    metrics_t m;
    float v0[3];
    int n, k, i;
    fxp_counters_t t;

    memset(&m, 0, sizeof m);
    s4_acc_init(&ap);
    s4_acc_init(&av);
    s4_acc_init(&aa);
    PROOF(dut->ctx_size <= S4_DUT_CTX_MAX, "DUT context fits the static buffer");
    s4_params_default(&prm);
    s4_dataset_filter_init(sc, v0);
    n = s4_dataset_generate(sc, DUR_S, g_frames, g_truth, S4_MAX_FRAMES);

    fxp_stats_reset();
    s4f_init(&mf, &prm);
    s4f_set_velocity(&mf, v0);
    memset(&g_ctx, 0, sizeof g_ctx);
    dut->init(&g_ctx, cfg, &prm, v0);

    for (k = 0; k < n; k++) {
        sensor_frame_fixed_t ff;
        es_ekf_fixed_t ex;
        double qx[4], qf[4], pd, nrm = 0.0;
        s4f_step(&mf, &g_frames[k]);
        sensor_frame_to_fixed(&g_frames[k], &ff, cfg);
        dut->step(&g_ctx, &ff);
        dut->export_state(&g_ctx, &ex);

        for (i = 0; i < 4; i++) {
            qx[i] = fxp_to_double_q(ex.x_nom.q[i], cfg->q.frac);
            qf[i] = mf.ekf.x_nom.q[i];
            nrm += qx[i] * qx[i];
        }
        if (fabs(sqrt(nrm) - 1.0) > m.qdev) m.qdev = fabs(sqrt(nrm) - 1.0);
        pd = 0.0;
        for (i = 0; i < 3; i++) {
            double p = fxp_to_double_q(ex.x_nom.p[i], cfg->p.frac);
            if (fabs(p) > pd) pd = fabs(p);
        }
        if (pd * 1000.0 > m.drift_mm) m.drift_mm = pd * 1000.0;

        if ((k + 1) % 50 == 0) {
            for (i = 0; i < 3; i++) {
                double px = fxp_to_double_q(ex.x_nom.p[i], cfg->p.frac);
                double vx = fxp_to_double_q(ex.x_nom.v[i], cfg->v.frac);
                if (!isfinite(px) || !isfinite(vx)) m.nan_flag = 1;
                s4_acc_add(&ap, mf.ekf.x_nom.p[i], px);
                s4_acc_add(&av, mf.ekf.x_nom.v[i], vx);
            }
            s4_acc_add(&aa, 0.0, s4_quat_angle(qf, qx));
        }
    }
    {
        es_ekf_fixed_t ex;
        int pos = 1;
        dut->export_state(&g_ctx, &ex);
        for (i = 0; i < ERR_DIM; i++)
            if (ex.D[i] <= 0) pos = 0;
        m.d_positive = pos;
    }
    t = fxp_stats_total();
    m.sat = t.sat;
    m.ovf = t.ovf;
    m.pd = dut->pd_clamps ? dut->pd_clamps(&g_ctx) : 0;
    m.rmse_p = s4_acc_rmse(&ap);
    m.rmse_v = s4_acc_rmse(&av);
    m.rmse_a = s4_acc_rmse(&aa);
    return m;
}

/* ------------------------------------------------------------ acceptance */
static const s4_scenario_t SCEN[4] = {S4_SC_STATIC, S4_SC_CONST_VEL, S4_SC_DYNAMIC, S4_SC_HIGH_RANGE};

/* returns 1 if accepted; reason (<=96 chars) names the first violated criterion */
static int accept(s4_scenario_t sc, const metrics_t *m, char *reason)
{
    reason[0] = 0;
    if (m->nan_flag) { strcpy(reason, "NaN/Inf in state"); return 0; }
    if (m->sat > S4_ACC_MAX_SATURATIONS) { sprintf(reason, "%u saturations", m->sat); return 0; }
    if (m->pd > S4_ACC_MAX_PD_CLAMPS || !m->d_positive) { strcpy(reason, "D lost positivity"); return 0; }
    if (sc == S4_SC_STATIC) {
        if (m->drift_mm > S4_ACC_STATIC_DRIFT_MM) { sprintf(reason, "static drift %.3g mm", m->drift_mm); return 0; }
        return 1;
    }
    if (m->rmse_p > S4_ACC_POS_RMSE_M) { sprintf(reason, "position RMSE %.2e m", m->rmse_p); return 0; }
    if (m->rmse_v > S4_ACC_VEL_RMSE_MPS) { sprintf(reason, "velocity RMSE %.2e m/s", m->rmse_v); return 0; }
    if (m->rmse_a > S4_ACC_ATT_RMSE_RAD) { sprintf(reason, "attitude RMSE %.2e rad", m->rmse_a); return 0; }
    if (m->qdev > S4_ACC_QNORM_DEV) { sprintf(reason, "|q| deviation %.2e", m->qdev); return 0; }
    return 1;
}

/* ------------------------------------------------------------ candidates */
typedef struct { const char *name; int dx, U, D, Phi; } cand_t;
static const cand_t CAND[] = {{"W16 (16/16/16/16)", 16, 16, 16, 16},
                              {"mixed dx16 U24 D32 Phi16", 16, 24, 32, 16},
                              {"W24 (24/24/24/24)", 24, 24, 24, 24},
                              {"W24 + D32 (dx24 U24 D32 Phi24)", 24, 24, 32, 24},
                              {"W24 + D32 + U32", 24, 32, 32, 24},
                              {"W32 (32/32/32/32)", 32, 32, 32, 32}};
#define NCAND ((int)(sizeof CAND / sizeof CAND[0]))

static fxp_fmt_t pick(int wl, int which)
{
    const ekf_qcfg_t *c = s4_cfg_get(wl == 16 ? 0 : (wl == 24 ? 1 : 2));
    return which == 0 ? c->dx : which == 1 ? c->U : which == 2 ? c->D : c->Phi;
}

static ekf_qcfg_t make_cfg(const cand_t *c)
{
    ekf_qcfg_t q = *s4_cfg_get(2); /* nominal path identical in every candidate */
    q.name = c->name;
    q.dx = pick(c->dx, 0);
    q.U = pick(c->U, 1);
    q.D = pick(c->D, 2);
    q.Phi = pick(c->Phi, 3);
    return q;
}

static int same_fmt(fxp_fmt_t a, fxp_fmt_t b) { return a.wl == b.wl && a.frac == b.frac; }

/* ------------------------------------------------------------------- main */
int main(int argc, char **argv)
{
    const char *dir = (argc > 1) ? argv[1] : "results";
    char path[1024];
    FILE *md, *csv;
    metrics_t M[NCAND][4];
    int ok[NCAND], order[NCAND], i, j, s, sel = -1;
    char why[NCAND][128];
    const ekf_qcfg_t *best = s4_best_cfg();

    printf("== System 4 proof of the best configuration (60 s synthetic scenarios)\n");

    /* ---------- A. dummy arithmetic checks ---------- */
    printf("-- A. dummy arithmetic checks\n");
    fxp_stats_reset();
    PROOF(fxp_add_sat(INT32_MAX, 1) == INT32_MAX && fxp_sub_sat(INT32_MIN, 1) == INT32_MIN, "add/sub saturate, no wrap-around");
    PROOF(fxp_mul(3 * 65536 / 2, 2 * 65536) == 3 * 65536, "Q16.16 1.5*2 = 3");
    PROOF(fxp_mul_round(-32768, 1) == -1 && fxp_mul(-32768, 1) == -1 && fxp_mul_round(32768, 1) == 1, "rounding symmetric, floor biased");
    PROOF(fxp_from_double_q(1.0, 16, 15) == 32767, "+1.0 saturates in Q1.15");
    PROOF(fxp_from_double_q(1.0, 32, 30) == (1 << 30), "+1.0 exact in Q2.30");
    PROOF(fxp_div(65536, 0) == INT32_MAX && fxp_stats_total().div0 == 1, "division by zero saturates and is counted");
    {
        int32_t sn, cs;
        double e = 0;
        for (i = 0; i < 2001; i++) {
            double a = -3.1 + 6.2 * i / 2000.0, aq = fxp_to_double_q(cordic_angle_from_double(a), 29);
            cordic_sincos(cordic_angle_from_double(a), &sn, &cs);
            e = fmax(e, fmax(fabs(fxp_to_double_q(sn, 30) - sin(aq)), fabs(fxp_to_double_q(cs, 30) - cos(aq))));
        }
        PROOF(e < 4e-8, "CORDIC error < 4e-8 against libm");
    }
    {
        int bad = 0;
        for (i = 0; i < 20000; i++) {
            int32_t a = (int32_t)((unsigned)rand() % 4000000u) - 2000000, b = (int32_t)((unsigned)rand() % 400000u) + 1;
            int32_t r0 = div_q(a, b, 16, 16, 16, DIV_BASELINE, NULL);
            if (div_q(a, b, 16, 16, 16, DIV_RESTORING, NULL) != r0 || div_q(a, b, 16, 16, 16, DIV_NONRESTORING, NULL) != r0) bad++;
        }
        PROOF(bad == 0, "restoring and non-restoring dividers bit-exact with baseline");
    }

    /* ---------- B. candidate search ---------- */
    printf("-- B. candidate search (storage cost ascending)\n");
    for (i = 0; i < NCAND; i++) {
        ekf_qcfg_t c = make_cfg(&CAND[i]);
        ok[i] = 1;
        why[i][0] = 0;
        for (s = 0; s < 4; s++) {
            char r[96];
            M[i][s] = run_dut(&S4_DUT_STUDY, &c, SCEN[s]);
            if (!accept(SCEN[s], &M[i][s], r) && ok[i]) {
                ok[i] = 0;
                snprintf(why[i], sizeof why[i], "%s: %s", s4_scenario_name(SCEN[s]), r);
            }
        }
        order[i] = i;
    }
    for (i = 0; i < NCAND; i++) /* sort by storage bits (insertion sort, fixed small N) */
        for (j = i + 1; j < NCAND; j++) {
            ekf_qcfg_t a = make_cfg(&CAND[order[i]]), b = make_cfg(&CAND[order[j]]);
            if (s4_cfg_storage_bits(&b) < s4_cfg_storage_bits(&a)) { int t = order[i]; order[i] = order[j]; order[j] = t; }
        }
    for (i = 0; i < NCAND; i++) {
        ekf_qcfg_t c = make_cfg(&CAND[order[i]]);
        printf("   %-34s %6u bit  %s %s\n", CAND[order[i]].name, s4_cfg_storage_bits(&c),
               ok[order[i]] ? "ACCEPTED" : "rejected", ok[order[i]] ? "" : why[order[i]]);
        if (sel < 0 && ok[order[i]]) sel = order[i];
    }
    PROOF(sel >= 0, "at least one candidate meets every acceptance threshold");
    PROOF(!ok[0] || sel == order[0], "selection is the cheapest accepted candidate");

    /* ---------- C. header consistency ---------- */
    printf("-- C. header consistency\n");
    if (sel >= 0) {
        ekf_qcfg_t c = make_cfg(&CAND[sel]);
        printf("   selected: %s | s4_best_cfg(): dx Q%d.%d U Q%d.%d D Q%d.%d Phi Q%d.%d\n", CAND[sel].name,
               best->dx.wl - best->dx.frac, best->dx.frac, best->U.wl - best->U.frac, best->U.frac,
               best->D.wl - best->D.frac, best->D.frac, best->Phi.wl - best->Phi.frac, best->Phi.frac);
        PROOF(same_fmt(c.dx, best->dx) && same_fmt(c.U, best->U) && same_fmt(c.D, best->D) && same_fmt(c.Phi, best->Phi),
              "selected configuration equals s4_best_cfg() in system4.h");
        PROOF(same_fmt(c.p, best->p) && same_fmt(c.v, best->v) && same_fmt(c.q, best->q) && same_fmt(c.b, best->b),
              "nominal-path formats equal s4_best_cfg()");
    }
    {
        metrics_t mb;
        char r[96];
        int all = 1;
        for (s = 0; s < 4; s++) {
            mb = run_dut(&S4_DUT_STUDY, best, SCEN[s]);
            if (!accept(SCEN[s], &mb, r)) all = 0;
        }
        PROOF(all, "s4_best_cfg() itself passes all four scenarios");
    }

    /* ---------- D. negative control ---------- */
    printf("-- D. negative control: TODO template must be rejected\n");
    {
        int rejected = 0;
        for (s = 0; s < 4; s++) {
            metrics_t mt = run_dut(&S4_DUT_TEMPLATE, best, SCEN[s]);
            char r[96];
            int a = accept(SCEN[s], &mt, r);
            printf("   template on %-9s %s %s\n", s4_scenario_name(SCEN[s]), a ? "accepted" : "REJECTED", a ? "" : r);
            if (!a) rejected = 1;
        }
        PROOF(rejected, "harness rejects the unimplemented template (it can fail)");
    }

    /* ---------- E. report ---------- */
    snprintf(path, sizeof path, "%s/proof/best_config_proof.md", dir);
    md = fopen(path, "w");
    snprintf(path, sizeof path, "%s/proof/best_config_proof.csv", dir);
    csv = fopen(path, "w");
    if (!md || !csv) {
        fprintf(stderr, "cannot write report under %s/proof (create the directory)\n", dir);
        return 2;
    }
    fprintf(csv, "candidate,storage_bits,scenario,rmse_pos_m,rmse_vel_mps,rmse_att_rad,qnorm_dev,static_drift_mm,saturations,accepted\n");
    fprintf(md, "# System 4 -- best-configuration proof (generated by tests/proof/system4_proof.c)\n\n"
                "Synthetic data, 60 s, fixed vs float32 reference, study pipeline (PENDING INTEGRATION).\n"
                "Thresholds (PROPOSED, include/system4.h): position %.0e m, velocity %.0e m/s, attitude %.0e rad, "
                "|q| deviation %.0e, saturations %d, static drift %.1f mm, D>0.\n\n"
                "| candidate | storage bits | verdict | first violation |\n|---|---|---|---|\n",
            S4_ACC_POS_RMSE_M, S4_ACC_VEL_RMSE_MPS, S4_ACC_ATT_RMSE_RAD, S4_ACC_QNORM_DEV, S4_ACC_MAX_SATURATIONS,
            S4_ACC_STATIC_DRIFT_MM);
    for (i = 0; i < NCAND; i++) {
        ekf_qcfg_t c = make_cfg(&CAND[order[i]]);
        fprintf(md, "| %s%s | %u | %s | %s |\n", CAND[order[i]].name, order[i] == sel ? " **(selected)**" : "",
                s4_cfg_storage_bits(&c), ok[order[i]] ? "ACCEPTED" : "rejected", why[order[i]]);
        for (s = 0; s < 4; s++) {
            char r[96];
            const metrics_t *m = &M[order[i]][s];
            fprintf(csv, "%s,%u,%s,%.6e,%.6e,%.6e,%.6e,%.6f,%u,%d\n", CAND[order[i]].name, s4_cfg_storage_bits(&c),
                    s4_scenario_name(SCEN[s]), m->rmse_p, m->rmse_v, m->rmse_a, m->qdev, m->drift_mm, m->sat,
                    accept(SCEN[s], m, r));
        }
    }
    if (sel >= 0) {
        fprintf(md, "\n## Selected configuration, per scenario\n\n| scenario | pos RMSE [m] | vel RMSE [m/s] | att RMSE [rad] | max abs(|q|-1) | saturations |\n|---|---|---|---|---|---|\n");
        for (s = 0; s < 4; s++)
            fprintf(md, "| %s | %.2e | %.2e | %.2e | %.2e | %u |\n", s4_scenario_name(SCEN[s]), M[sel][s].rmse_p,
                    M[sel][s].rmse_v, M[sel][s].rmse_a, M[sel][s].qdev, M[sel][s].sat);
    }
    fprintf(md, "\nChecks: %d, failed: %d.\n", g_checks, g_fails);
    fclose(md);
    fclose(csv);

    printf("== %d proof checks, %d failed  ->  %s/proof/best_config_proof.md\n", g_checks, g_fails, dir);
    return g_fails ? 1 : 0;
}
