/*
 * system4_dut.c -- DUT bindings declared in system4.h.
 *   S4_DUT_STUDY    : wraps the study pipeline (ekf_s4_ref_fixed.c)
 *   S4_DUT_TEMPLATE : algorithm template with TODO stages (see TEMPLATE markers).
 *                     Copy this block, fill each stage, and run build/bin/system4_proof.
 */
#include <math.h> /* init-time conversion only */
#include <string.h>

#include "system4.h"

/* ---------------------------------------------------------- study binding */
static void study_init(void *c, const ekf_qcfg_t *cfg, const s4_params_t *p, const float v0[3])
{
    s4x_init((s4x_model_t *)c, cfg, p);
    s4x_set_velocity((s4x_model_t *)c, v0);
}
static void study_step(void *c, const sensor_frame_fixed_t *f) { s4x_step((s4x_model_t *)c, f); }
static void study_export(const void *c, es_ekf_fixed_t *o) { *o = ((const s4x_model_t *)c)->ekf; }
static uint32_t study_pd(const void *c) { return ((const s4x_model_t *)c)->n_pd_clamp; }

const s4_dut_t S4_DUT_STUDY = {"study pipeline", sizeof(s4x_model_t), study_init, study_step,
                               study_export, study_pd};

/* ------------------------------------------------------ algorithm template */
typedef struct {
    es_ekf_fixed_t ekf; /* nominal (cfg formats) + normalized dx/U/D */
    ekf_qcfg_t cfg;
    int32_t dt_acc;     /* accumulated dt for the covariance step, Q2.30 */
    int cnt, decim;
} tpl_ctx_t;

/* TEMPLATE stage 1 -- nominal propagation (member 1).
 * In : x_nom, a_m/w_m (Q8.24), dt (Q2.30).  Out: x_nom updated, q normalised.
 * TODO(member 1): a_unb = a_m - ba ; an = R(q) a_unb + g ; p,v,q update (ICD Role 1). */
static void tpl_predict_nominal(tpl_ctx_t *t, const sensor_frame_fixed_t *f)
{
    (void)t;
    (void)f;
}

/* TEMPLATE stage 2 -- covariance time update (members 1+3).
 * In : U, D (normalized), dPhi (Q1.x, identity NOT stored), Qd.  Out: U, D.
 * TODO(member 3): Thornton modified weighted Gram-Schmidt, fixed loop bounds. */
static void tpl_covariance_update(tpl_ctx_t *t)
{
    (void)t;
}

/* TEMPLATE stage 3 -- one scalar measurement channel (members 2+3).
 * In : channel index i0 (H = e_i0), normalized innovation y_n (Q8.24), R_n (Q3.29).
 * TODO(member 3): Bierman scalar update of U, D and accumulation into dx. */
static void tpl_scalar_update(tpl_ctx_t *t, int i0, int32_t y_n, int32_t r_n)
{
    (void)t;
    (void)i0;
    (void)y_n;
    (void)r_n;
}

/* TEMPLATE stage 4 -- error injection and reset (member 2).
 * TODO(member 2): x_nom += dx ; q = normalize(q (x) [1, dtheta/2]) ; dx = 0. */
static void tpl_inject_reset(tpl_ctx_t *t)
{
    (void)t;
}

static void tpl_init(void *c, const ekf_qcfg_t *cfg, const s4_params_t *p, const float v0[3])
{
    tpl_ctx_t *t = (tpl_ctx_t *)c;
    int i, j;
    memset(t, 0, sizeof *t);
    t->cfg = *cfg;
    t->decim = p->cov_decim;
    t->ekf.x_nom.q[0] = (int32_t)((int64_t)1 << 30);
    for (i = 0; i < 3; i++) t->ekf.x_nom.v[i] = fxp_from_double_q((double)v0[i], cfg->v.wl, cfg->v.frac);
    for (j = 0; j < ERR_DIM; j++) {
        for (i = 0; i <= j; i++)
            t->ekf.U[U_IDX(i, j)] = (i == j) ? (int32_t)((int64_t)1 << cfg->U.frac) : 0;
        t->ekf.D[j] = fxp_from_double_q(ldexp((double)p->D0[j / 3], -2 * s4_scale_log2[j]), cfg->D.wl, cfg->D.frac);
    }
}

static void tpl_step(void *c, const sensor_frame_fixed_t *f)
{
    tpl_ctx_t *t = (tpl_ctx_t *)c;
    tpl_predict_nominal(t, f);
    t->dt_acc = fxp_add_sat(t->dt_acc, f->dt);
    if (++t->cnt >= t->decim) {
        tpl_covariance_update(t);
        t->dt_acc = 0;
        t->cnt = 0;
    }
    if (f->type == SENS_GNSS_POS_VEL || f->type == SENS_BARO) {
        tpl_covariance_update(t); /* flush pending dt before the measurement */
        tpl_scalar_update(t, 0, 0, 0); /* TODO: loop over the channels of this sensor */
        tpl_inject_reset(t);
    }
}

static void tpl_export(const void *c, es_ekf_fixed_t *o) { *o = ((const tpl_ctx_t *)c)->ekf; }

const s4_dut_t S4_DUT_TEMPLATE = {"algorithm template (TODO stubs)", sizeof(tpl_ctx_t), tpl_init,
                                  tpl_step, tpl_export, NULL};
