/*
 * ekf_s4_ref_fixed.c -- Model B: bit-accurate fixed-point study pipeline.
 *
 * Mirrors ekf_s4_ref_float.c operation by operation, in the NORMALIZED domain
 * for dx/U/D (see ekf_fixed_types.h).  Only integer arithmetic is used in the
 * core; the single block between BEGIN-CONVERSION / END-CONVERSION converts
 * the PROPOSED float parameters once at initialisation.
 *
 * Storage vs. compute: dx/U/D/Phi are STORED in the configured word length
 * (16/24/32 bit).  Intermediates use 32-bit working registers (Q8.24 W,
 * Q3.29 alpha/g/k) and 64-bit accumulators.  Whether a real datapath can
 * afford those widths is a DSP48/BRAM decision that is PENDING.
 *
 * Divisions (per (i,j) pair, not a shared reciprocal) go through the divider
 * module (div_set_pipeline_mode selects baseline/restoring/non-restoring;
 * all are bit-exact with each other).
 */
#include <math.h> /* BEGIN/END-CONVERSION block only */
#include <string.h>

#include "ekf_fixed_adapter.h"
#include "ekf_fixed_divider.h"
#include "ekf_fixed_math.h"
#include "ekf_s4_models.h"

#define QONE ((int64_t)1 << 30) /* 1.0 in Q2.30 */
#define WFRAC 24                /* S4_WORK_FRAC */
#define AFRAC 29                /* S4_ALPHA_FRAC */

static int32_t sat32(int64_t v) { return fxp_saturate(v); }

static int64_t div_pipe(int64_t n, int64_t d)
{
    div_info_t info;
    return divs64_round(n, d, div_get_pipeline_mode(), &info);
}

/* ------------------------------------------------------------ init */
/* BEGIN-CONVERSION */
void s4x_init(s4x_model_t *m, const ekf_qcfg_t *cfg, const s4_params_t *prm)
{
    int i, j, g;
    memset(m, 0, sizeof *m);
    m->cfg = *cfg;
    m->decim = prm->cov_decim;
    for (j = 0; j < ERR_DIM; j++) {
        for (i = 0; i <= j; i++)
            m->ekf.U[U_IDX(i, j)] =
                (i == j) ? fxp_sat_wl((int64_t)1 << cfg->U.frac, cfg->U.wl) : 0;
        m->ekf.D[j] = fxp_from_double_q(
            ldexp((double)prm->D0[j / 3], -2 * s4_scale_log2[j]), cfg->D.wl, cfg->D.frac);
    }
    m->ekf.x_nom.q[0] = (int32_t)QONE;
    for (g = 0; g < 5; g++) {
        double v = ldexp((double)prm->qc[g], -2 * s4_scale_log2[3 * g]);
        m->qc_n[g] = (int64_t)(ldexp(v, 50) + 0.5);
    }
    m->r_n[0] = fxp_from_double_q(ldexp((double)prm->r_pos, -2 * s4_scale_log2[0]), 32u, AFRAC);
    m->r_n[1] = fxp_from_double_q(ldexp((double)prm->r_vel, -2 * s4_scale_log2[3]), 32u, AFRAC);
    m->r_n[2] = fxp_from_double_q(ldexp((double)prm->r_baro, -2 * s4_scale_log2[2]), 32u, AFRAC);
    m->g_q = fxp_from_double_q(S4_GRAVITY_MPS2, 32u, S4_IMU_FRAC);
    m->Rm[0] = m->Rm[4] = m->Rm[8] = (int32_t)QONE;
}

void s4x_set_velocity(s4x_model_t *m, const float v[3])
{
    int i;
    for (i = 0; i < 3; i++)
        m->ekf.x_nom.v[i] = fxp_from_double_q((double)v[i], m->cfg.v.wl, m->cfg.v.frac);
}
/* END-CONVERSION */

/* ------------------------------------------------- quaternion helpers */
static void rot_from_q(const int32_t *q, int32_t *R)
{
    int64_t w = q[0], x = q[1], y = q[2], z = q[3];
#define PQ(a, b) fxp_rshift_round((a) * (b), 30)
    R[0] = sat32(QONE - 2 * (PQ(y, y) + PQ(z, z)));
    R[1] = sat32(2 * (PQ(x, y) - PQ(w, z)));
    R[2] = sat32(2 * (PQ(x, z) + PQ(w, y)));
    R[3] = sat32(2 * (PQ(x, y) + PQ(w, z)));
    R[4] = sat32(QONE - 2 * (PQ(x, x) + PQ(z, z)));
    R[5] = sat32(2 * (PQ(y, z) - PQ(w, x)));
    R[6] = sat32(2 * (PQ(x, z) - PQ(w, y)));
    R[7] = sat32(2 * (PQ(y, z) + PQ(w, x)));
    R[8] = sat32(QONE - 2 * (PQ(x, x) + PQ(y, y)));
#undef PQ
}

static void qmul_fx(const int32_t *p, const int32_t *q, int32_t *r)
{
    int64_t p0 = p[0], p1 = p[1], p2 = p[2], p3 = p[3];
    int64_t q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
    r[0] = sat32(fxp_rshift_round(p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3, 30));
    r[1] = sat32(fxp_rshift_round(p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2, 30));
    r[2] = sat32(fxp_rshift_round(p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1, 30));
    r[3] = sat32(fxp_rshift_round(p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0, 30));
}

/* Quaternion normalisation.
 * Default: inverse-sqrt by Newton-Raphson, y <- y (3 - n2 y^2) / 2, started from a
 * power-of-two estimate (leading-one position) so n2 y0^2 is in [0.5, 2) and 4 fixed
 * iterations converge; no sqrt, no division.  Compile with -DS4_QNORM_TAYLOR for the
 * ICD order-1 rule q *= (3 - |q|^2)/2, which is only valid for |q| very close to 1. */
static void qnorm_fx(int32_t *q)
{
    int64_t n2 = 0, fac;
    int i;
    for (i = 0; i < 4; i++) n2 += (int64_t)q[i] * q[i];
    n2 = fxp_rshift_round(n2, 30); /* Q2.30 */
#ifdef S4_QNORM_TAYLOR
    fac = ((int64_t)3 * QONE - n2) >> 1;
#else
    if (n2 <= 0) return;
    {
        int p = 0, e, it;
        int64_t y;
        while ((n2 >> (p + 1)) != 0 && p < 62) p++; /* index of leading one */
        e = p - 30;                                  /* floor(log2 value) */
        y = (e + 1 >= 0) ? (QONE >> ((e + 1) >> 1)) : (QONE << (-((e + 1) >> 1)));
        for (it = 0; it < 4; it++) {
            int64_t t = fxp_rshift_round(y * y, 30);
            int64_t u = fxp_rshift_round(n2 * t, 30);
            y = fxp_rshift_round(y * ((int64_t)3 * QONE - u), 31);
        }
        fac = y;
    }
#endif
    for (i = 0; i < 4; i++) q[i] = sat32(fxp_rshift_round((int64_t)q[i] * fac, 30));
}

/* -------------------------------------------------- nominal propagation */
static void nominal_step(s4x_model_t *m, const sensor_frame_fixed_t *f)
{
    const ekf_qcfg_t *c = &m->cfg;
    nominal_state_fixed_t *x = &m->ekf.x_nom;
    int32_t an[3], adt[3], dq[4], qn[4];
    int32_t dt = f->dt;
    int i, j;

    fxp_stats_set_group(FXP_GRP_NOM_B);
    for (i = 0; i < 3; i++) {
        int64_t ba24 = fxp_rshift_round(x->ba[i], (int)c->b.frac - (int)S4_IMU_FRAC);
        int64_t bg24 = fxp_rshift_round(x->bg[i], (int)c->b.frac - (int)S4_IMU_FRAC);
        m->a_unb[i] = sat32((int64_t)f->imu_acc[i] - ba24);
        m->w_unb[i] = sat32((int64_t)f->imu_gyr[i] - bg24);
    }
    fxp_stats_set_group(FXP_GRP_NOM_Q);
    rot_from_q(x->q, m->Rm);
    fxp_stats_set_group(FXP_GRP_NOM_V);
    for (i = 0; i < 3; i++) {
        int64_t s = 0;
        for (j = 0; j < 3; j++) s += (int64_t)m->Rm[3 * i + j] * m->a_unb[j]; /* Q54 */
        an[i] = sat32(fxp_rshift_round(s, 30));                               /* Q8.24 */
    }
    an[2] = sat32((int64_t)an[2] + m->g_q);
    for (i = 0; i < 3; i++) adt[i] = sat32(fxp_rshift_round((int64_t)an[i] * dt, 30));

    for (i = 0; i < 3; i++) {
        int64_t t1, t2, inc;
        fxp_stats_set_group(FXP_GRP_NOM_P);
        t1 = fxp_rshift_round((int64_t)x->v[i] * dt, (int)c->v.frac + 30 - 54); /* frac 54 */
        t2 = fxp_asr64((int64_t)adt[i] * dt, 1);
        inc = fxp_rshift_round(t1 + t2, 54 - (int)c->p.frac);
        if (inc == 0 && (t1 + t2) != 0) fxp_stats_note_underflow(); /* increment lost below 1 LSB */
        x->p[i] = fxp_sat_wl((int64_t)x->p[i] + inc, c->p.wl);
        fxp_stats_set_group(FXP_GRP_NOM_V);
        x->v[i] = fxp_sat_wl((int64_t)x->v[i] +
                                 fxp_rshift_round(adt[i], (int)S4_IMU_FRAC - (int)c->v.frac),
                             c->v.wl);
    }

    fxp_stats_set_group(FXP_GRP_NOM_Q);
    dq[0] = (int32_t)QONE;
    for (i = 0; i < 3; i++)
        dq[1 + i] = sat32(fxp_rshift_round((int64_t)m->w_unb[i] * dt, 25)); /* 1/2 w dt, Q2.30 */
    qmul_fx(x->q, dq, qn);
    qnorm_fx(qn);
    memcpy(x->q, qn, sizeof qn);
}

/* ------------------------------------------------------------ Phi (dPhi) */
/* Normalised dPhi = (Phi - I) scaled: dPhi_n(i,j) = F(i,j)*dt*2^(e_j-e_i).
 * The identity part is NOT stored: W = U + dPhi*U is formed exactly.      */
static void build_phi_fx(const s4x_model_t *m, int32_t dtc, int32_t *phi)
{
    const ekf_qcfg_t *c = &m->cfg;
    int fp = (int)c->Phi.frac;
    uint8_t wl = c->Phi.wl;
    int64_t sa[9], sw[9];
    const int32_t *a = m->a_unb, *w = m->w_unb;
    int i, j, k;

    fxp_stats_set_group(FXP_GRP_PHI);
    memset(phi, 0, sizeof(int32_t) * ERR_DIM * ERR_DIM);
    sa[0] = 0; sa[1] = -(int64_t)a[2]; sa[2] = a[1];
    sa[3] = a[2]; sa[4] = 0; sa[5] = -(int64_t)a[0];
    sa[6] = -(int64_t)a[1]; sa[7] = a[0]; sa[8] = 0;
    sw[0] = 0; sw[1] = -(int64_t)w[2]; sw[2] = w[1];
    sw[3] = w[2]; sw[4] = 0; sw[5] = -(int64_t)w[0];
    sw[6] = -(int64_t)w[1]; sw[7] = w[0]; sw[8] = 0;

    for (i = 0; i < 3; i++) {
        phi[MAT15_IDX(i, 3 + i)] = fxp_sat_wl(fxp_rshift_round(dtc, 30 - fp), wl);
        phi[MAT15_IDX(6 + i, 12 + i)] = fxp_sat_wl(-fxp_rshift_round(dtc, 30 - fp + 1), wl);
        for (j = 0; j < 3; j++) {
            int64_t ra = 0, m24;
            for (k = 0; k < 3; k++) ra += (int64_t)m->Rm[3 * i + k] * sa[3 * k + j]; /* Q54 */
            m24 = fxp_rshift_round(ra, 30);                                         /* Q8.24 */
            phi[MAT15_IDX(3 + i, 6 + j)] =
                fxp_sat_wl(fxp_rshift_round(-m24 * (int64_t)dtc, 54 - fp + 5), wl);
            phi[MAT15_IDX(3 + i, 9 + j)] = fxp_sat_wl(
                fxp_rshift_round(-(int64_t)m->Rm[3 * i + j] * (int64_t)dtc, 60 - fp + 3), wl);
            phi[MAT15_IDX(6 + i, 6 + j)] =
                fxp_sat_wl(fxp_rshift_round(-sw[3 * i + j] * (int64_t)dtc, 54 - fp), wl);
        }
    }
}

/* ----------------------------------------------------- Thornton (fixed) */
/* sum_k D_k W_ik W_jk + sum_k Qd_k W_i,15+k W_j,15+k ; result frac = fd+WFRAC */
static int64_t dotD(const int32_t *wi, const int32_t *wj, const int32_t *D, const int64_t *qd)
{
    int64_t acc = 0;
    int k;
    for (k = 0; k < ERR_DIM; k++) {
        if (wi[k] != 0 && wj[k] != 0) {
            int64_t t = fxp_rshift_round((int64_t)wi[k] * wj[k], WFRAC);
            acc = fxp_add64_sat(acc, fxp_mul64_sat((int64_t)D[k], t));
        }
    }
    for (k = 0; k < ERR_DIM; k++) {
        if (wi[ERR_DIM + k] != 0 && wj[ERR_DIM + k] != 0) {
            int64_t t = fxp_rshift_round((int64_t)wi[ERR_DIM + k] * wj[ERR_DIM + k], WFRAC);
            acc = fxp_add64_sat(acc, fxp_rshift_round(fxp_mul64_sat(qd[k], t), WFRAC));
        }
    }
    return acc;
}

/* value num/den (den > 0) as raw Qfr, saturated to wl bits */
static int32_t ratio_to_q(int64_t num, int64_t den, int fr, uint8_t wl)
{
    int sh = 0;
    int64_t n, d, lim, q;
    while ((den >> sh) >= ((int64_t)1 << 31)) sh++; /* <= 32 iterations */
    d = den >> sh;
    n = fxp_asr64(num, sh);
    lim = INT64_MAX >> fr;
    if (n > lim) return fxp_sat_wl(INT64_MAX, wl);
    if (n < -lim) return fxp_sat_wl(INT64_MIN, wl);
    q = div_pipe(n * ((int64_t)1 << fr), d);
    return fxp_sat_wl(q, wl);
}

static int32_t store_D(s4x_model_t *m, int64_t raw)
{
    fxp_stats_set_group(FXP_GRP_D);
    if (raw <= 0) {
        m->n_pd_clamp++;
        raw = 1;
    }
    return fxp_sat_wl(raw, m->cfg.D.wl);
}

static void thornton_fx(s4x_model_t *m, const int32_t *phi, const int64_t *qd)
{
    const ekf_qcfg_t *c = &m->cfg;
    int32_t *U = m->ekf.U;
    int32_t *D = m->ekf.D;
    int fu = (int)c->U.frac, fp = (int)c->Phi.frac;
    int32_t W[ERR_DIM][2 * ERR_DIM];
    int32_t Dn[ERR_DIM];
    int32_t Un[U_UPPER_SIZE];
    int i, j, k;

    fxp_stats_set_group(FXP_GRP_WORK);
    for (i = 0; i < ERR_DIM; i++) {
        for (j = 0; j < ERR_DIM; j++) {
            int64_t acc = 0, w;
            for (k = 0; k <= j; k++) {
                int64_t ph = phi[MAT15_IDX(i, k)];
                if (ph != 0) acc = fxp_add64_sat(acc, ph * (int64_t)U[U_IDX(k, j)]);
            }
            w = fxp_rshift_round(acc, fp + fu - WFRAC);
            if (i <= j) w += fxp_rshift_round(U[U_IDX(i, j)], fu - WFRAC);
            W[i][j] = sat32(w);
            W[i][ERR_DIM + j] = (i == j) ? (int32_t)((int64_t)1 << WFRAC) : 0;
        }
    }

    for (j = ERR_DIM - 1; j >= 0; j--) {
        int64_t accD = dotD(W[j], W[j], D, qd);
        Dn[j] = store_D(m, fxp_rshift_round(accD, WFRAC));
        Un[U_IDX(j, j)] = fxp_sat_wl((int64_t)1 << fu, c->U.wl);
        for (i = 0; i < j; i++) {
            int32_t u;
            fxp_stats_set_group(FXP_GRP_U);
            if (accD <= 0) {
                u = 0;
            } else {
                u = ratio_to_q(dotD(W[i], W[j], D, qd), accD, fu, c->U.wl);
            }
            Un[U_IDX(i, j)] = u;
            fxp_stats_set_group(FXP_GRP_WORK);
            for (k = 0; k < 2 * ERR_DIM; k++) {
                if (W[j][k] != 0)
                    W[i][k] = sat32((int64_t)W[i][k] -
                                    fxp_rshift_round((int64_t)u * W[j][k], fu));
            }
        }
    }
    memcpy(U, Un, sizeof Un);
    memcpy(D, Dn, sizeof Dn);
}

static void cov_update_fx(s4x_model_t *m)
{
    int32_t phi[ERR_DIM * ERR_DIM];
    int64_t qd[ERR_DIM];
    int fd = (int)m->cfg.D.frac;
    int i;
    if (m->dt_acc <= 0) return;
    build_phi_fx(m, m->dt_acc, phi);
    for (i = 0; i < ERR_DIM; i++) /* qc_n frac 50 * dt frac 30 = 80 -> fd+WFRAC */
        qd[i] = fxp_rshift_round(fxp_mul64_sat(m->qc_n[i / 3], (int64_t)m->dt_acc), 80 - (fd + WFRAC));
    thornton_fx(m, phi, qd);
    m->dt_acc = 0;
    m->cnt = 0;
}

/* ------------------------------------------------------ Bierman (fixed) */
/* y_n : normalised innovation Q8.24 ; R_n : Q3.29 */
static void bierman_fx(s4x_model_t *m, int i0, int32_t y_n, int32_t R_n)
{
    const ekf_qcfg_t *c = &m->cfg;
    int32_t *U = m->ekf.U;
    int32_t *D = m->ekf.D;
    int32_t *dx = m->ekf.dx;
    int fu = (int)c->U.frac, fd = (int)c->D.frac, fdx = (int)c->dx.frac;
    int32_t f[ERR_DIM], g[ERR_DIM], kv[ERR_DIM];
    int32_t alpha = R_n, innov;
    int i, j;

    fxp_stats_set_group(FXP_GRP_WORK);
    for (j = 0; j < ERR_DIM; j++) {
        kv[j] = 0;
        f[j] = (j >= i0) ? U[U_IDX(i0, j)] : 0;
        g[j] = sat32(fxp_rshift_round((int64_t)D[j] * f[j], fd + fu - AFRAC));
    }
    for (j = i0; j < ERR_DIM; j++) {
        int32_t an = sat32((int64_t)alpha + fxp_rshift_round((int64_t)f[j] * g[j], fu));
        int64_t dnew;
        if (an <= 0) { /* cannot happen with R>0 unless saturated */
            m->n_bierman_abort++;
            return;
        }
        dnew = div_pipe(fxp_mul64_sat((int64_t)D[j], (int64_t)alpha), an); /* frac fd */
        kv[j] = g[j];
        for (i = 0; i < j; i++) {
            int32_t tmp = U[U_IDX(i, j)];
            int64_t corr = div_pipe(fxp_mul64_sat((int64_t)f[j], (int64_t)kv[i]), alpha); /* frac fu */
            fxp_stats_set_group(FXP_GRP_U);
            U[U_IDX(i, j)] = fxp_sat_wl((int64_t)tmp - corr, c->U.wl);
            fxp_stats_set_group(FXP_GRP_WORK);
            kv[i] = sat32((int64_t)kv[i] + fxp_rshift_round((int64_t)tmp * g[j], fu));
        }
        D[j] = store_D(m, dnew);
        fxp_stats_set_group(FXP_GRP_WORK);
        alpha = an;
    }
    innov = sat32((int64_t)y_n - fxp_rshift_round(dx[i0], fdx - WFRAC));
    for (j = 0; j < ERR_DIM; j++) {
        int64_t num, delta;
        if (kv[j] == 0) continue;
        fxp_stats_set_group(FXP_GRP_WORK);
        num = fxp_mul64_sat((int64_t)kv[j], (int64_t)innov); /* frac 29+26 */
        num = fxp_rshift_round(num, WFRAC - fdx);
        delta = div_pipe(num, alpha);
        fxp_stats_set_group(FXP_GRP_DX);
        dx[j] = fxp_sat_wl((int64_t)dx[j] + delta, c->dx.wl);
    }
}

/* ---------------------------------------------------- inject and reset */
static void inject_reset_fx(s4x_model_t *m)
{
    const ekf_qcfg_t *c = &m->cfg;
    nominal_state_fixed_t *x = &m->ekf.x_nom;
    int32_t *dx = m->ekf.dx;
    int fdx = (int)c->dx.frac;
    int32_t dq[4], qn[4];
    int i;
    for (i = 0; i < 3; i++) {
        fxp_stats_set_group(FXP_GRP_NOM_P);
        x->p[i] = fxp_sat_wl((int64_t)x->p[i] + fxp_rshift_round(dx[i], fdx - ((int)c->p.frac + s4_scale_log2[i])), c->p.wl);
        fxp_stats_set_group(FXP_GRP_NOM_V);
        x->v[i] = fxp_sat_wl((int64_t)x->v[i] + fxp_rshift_round(dx[3 + i], fdx - ((int)c->v.frac + s4_scale_log2[3 + i])), c->v.wl);
        fxp_stats_set_group(FXP_GRP_NOM_B);
        x->ba[i] = fxp_sat_wl((int64_t)x->ba[i] + fxp_rshift_round(dx[9 + i], fdx - ((int)c->b.frac + s4_scale_log2[9 + i])), c->b.wl);
        x->bg[i] = fxp_sat_wl((int64_t)x->bg[i] + fxp_rshift_round(dx[12 + i], fdx - ((int)c->b.frac + s4_scale_log2[12 + i])), c->b.wl);
    }
    fxp_stats_set_group(FXP_GRP_NOM_Q);
    dq[0] = (int32_t)QONE;
    /* 1/2 * dtheta = dx_n * 2^(e-1) ; Q2.30 raw = dx_raw * 2^(30+e-1-fdx) */
    for (i = 0; i < 3; i++)
        dq[1 + i] = sat32(fxp_rshift_round(dx[6 + i], fdx - (30 + s4_scale_log2[6 + i] - 1)));
    qmul_fx(x->q, dq, qn);
    qnorm_fx(qn);
    memcpy(x->q, qn, sizeof qn);
    for (i = 0; i < ERR_DIM; i++) dx[i] = 0;
}

/* innovation (measurement - nominal) in the nominal format -> Q8.24 normalised */
static int32_t innovation_n(int64_t z, int64_t h, int frac, int state_idx)
{
    fxp_stats_set_group(FXP_GRP_WORK);
    return sat32(fxp_rshift_round(z - h, frac + s4_scale_log2[state_idx] - WFRAC));
}

/* ------------------------------------------------------------ step */
void s4x_step(s4x_model_t *m, const sensor_frame_fixed_t *f)
{
    const ekf_qcfg_t *c = &m->cfg;
    nominal_state_fixed_t *x = &m->ekf.x_nom;
    int i;

    nominal_step(m, f);

    fxp_stats_set_group(FXP_GRP_WORK);
    m->dt_acc = sat32((int64_t)m->dt_acc + f->dt);
    if (++m->cnt >= m->decim) cov_update_fx(m);

    if (f->type == SENS_GNSS_POS_VEL) {
        cov_update_fx(m);
        for (i = 0; i < 3; i++)
            bierman_fx(m, i, innovation_n(f->meas[i], x->p[i], (int)c->p.frac, i), m->r_n[0]);
        for (i = 0; i < 3; i++)
            bierman_fx(m, 3 + i, innovation_n(f->meas[3 + i], x->v[i], (int)c->v.frac, 3 + i),
                       m->r_n[1]);
        inject_reset_fx(m);
    } else if (f->type == SENS_BARO) {
        cov_update_fx(m);
        bierman_fx(m, 2, innovation_n(f->meas[0], x->p[2], (int)c->p.frac, 2), m->r_n[2]);
        inject_reset_fx(m);
    }
    /* SENS_MAG: PENDING */
}
