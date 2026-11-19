/*
 * ekf_s4_ref_float.c -- Model A: float32 golden reference of the study
 * pipeline.  Physical units, un-normalized, standard math (sqrtf for the
 * quaternion norm).  See ekf_s4_models.h for scope.
 */
#include <math.h>
#include <string.h>

#include "ekf_fixed_config.h"
#include "ekf_s4_models.h"

void s4_params_default(s4_params_t *p)
{
    static const float d0[5] = {1.0f, 0.5f, 0.05f, 0.01f, 0.0001f};     /* ICD */
    static const float qc[5] = {0.0f, 2.5e-3f, 2.5e-5f, 1.0e-6f, 1.0e-10f}; /* PROPOSED */
    int i;
    for (i = 0; i < 5; i++) {
        p->D0[i] = d0[i];
        p->qc[i] = qc[i];
    }
    p->r_pos = 0.25f;  /* ICD */
    p->r_vel = 0.04f;  /* ICD */
    p->r_baro = 0.16f; /* ICD */
    p->cov_decim = S4_COV_DECIM_DEFAULT;
}

static void quat_to_R(const float *q, float *R)
{
    float w = q[0], x = q[1], y = q[2], z = q[3];
    R[0] = 1.0f - 2.0f * (y * y + z * z);
    R[1] = 2.0f * (x * y - w * z);
    R[2] = 2.0f * (x * z + w * y);
    R[3] = 2.0f * (x * y + w * z);
    R[4] = 1.0f - 2.0f * (x * x + z * z);
    R[5] = 2.0f * (y * z - w * x);
    R[6] = 2.0f * (x * z - w * y);
    R[7] = 2.0f * (y * z + w * x);
    R[8] = 1.0f - 2.0f * (x * x + y * y);
}

static void quat_mul(const float *p, const float *q, float *r)
{
    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
    r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
    r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
    r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

static void quat_normalize(float *q)
{
    float n = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (n > 0.0f) {
        int i;
        for (i = 0; i < 4; i++) q[i] /= n;
    }
}

void s4f_init(s4f_model_t *m, const s4_params_t *prm)
{
    int i, j;
    memset(m, 0, sizeof *m);
    m->prm = *prm;
    m->ekf.x_nom.q[0] = 1.0f;
    for (j = 0; j < ERR_DIM; j++) {
        for (i = 0; i <= j; i++) m->ekf.U[U_IDX(i, j)] = (i == j) ? 1.0f : 0.0f;
        m->ekf.D[j] = prm->D0[j / 3];
    }
    m->R[0] = m->R[4] = m->R[8] = 1.0f;
}

void s4f_set_velocity(s4f_model_t *m, const float v[3])
{
    int i;
    for (i = 0; i < 3; i++) m->ekf.x_nom.v[i] = v[i];
}

/* Thornton modified weighted Gram-Schmidt: P+ = Phi U D U^T Phi^T + diag(qd) */
static uint32_t thornton(float *U, float *D, const float *Phi, const float *qd)
{
    float W[ERR_DIM][2 * ERR_DIM];
    float Dh[2 * ERR_DIM];
    float c[2 * ERR_DIM];
    float Dn[ERR_DIM];
    float Un[U_UPPER_SIZE];
    uint32_t viol = 0;
    int i, j, k;

    for (i = 0; i < ERR_DIM; i++) {
        for (j = 0; j < ERR_DIM; j++) {
            float acc = 0.0f;
            for (k = 0; k <= j; k++) acc += Phi[MAT15_IDX(i, k)] * U[U_IDX(k, j)];
            W[i][j] = acc;
            W[i][ERR_DIM + j] = (i == j) ? 1.0f : 0.0f;
        }
        Dh[i] = D[i];
        Dh[ERR_DIM + i] = qd[i];
    }
    for (j = ERR_DIM - 1; j >= 0; j--) {
        float dj = 0.0f;
        for (k = 0; k < 2 * ERR_DIM; k++) {
            c[k] = Dh[k] * W[j][k];
            dj += W[j][k] * c[k];
        }
        if (!(dj > 0.0f)) {
            viol++;
            dj = 1.0e-30f;
        }
        Dn[j] = dj;
        Un[U_IDX(j, j)] = 1.0f;
        for (i = 0; i < j; i++) {
            float s = 0.0f, u;
            for (k = 0; k < 2 * ERR_DIM; k++) s += W[i][k] * c[k];
            u = s / dj;
            Un[U_IDX(i, j)] = u;
            for (k = 0; k < 2 * ERR_DIM; k++) W[i][k] -= u * W[j][k];
        }
    }
    memcpy(U, Un, sizeof Un);
    memcpy(D, Dn, sizeof Dn);
    return viol;
}

static void build_phi(const s4f_model_t *m, float dt, float *Phi)
{
    const float *R = m->R, *a = m->a_unb, *w = m->w_unb;
    float Sa[9], Sw[9];
    int i, j, k;
    memset(Phi, 0, sizeof(float) * ERR_DIM * ERR_DIM);
    for (i = 0; i < ERR_DIM; i++) Phi[MAT15_IDX(i, i)] = 1.0f;
    Sa[0] = 0; Sa[1] = -a[2]; Sa[2] = a[1];
    Sa[3] = a[2]; Sa[4] = 0; Sa[5] = -a[0];
    Sa[6] = -a[1]; Sa[7] = a[0]; Sa[8] = 0;
    Sw[0] = 0; Sw[1] = -w[2]; Sw[2] = w[1];
    Sw[3] = w[2]; Sw[4] = 0; Sw[5] = -w[0];
    Sw[6] = -w[1]; Sw[7] = w[0]; Sw[8] = 0;
    for (i = 0; i < 3; i++) {
        Phi[MAT15_IDX(i, 3 + i)] = dt;             /* I dt           */
        Phi[MAT15_IDX(6 + i, 12 + i)] = -dt;       /* -I dt (gyro)   */
        for (j = 0; j < 3; j++) {
            float ra = 0.0f;
            for (k = 0; k < 3; k++) ra += R[3 * i + k] * Sa[3 * k + j];
            Phi[MAT15_IDX(3 + i, 6 + j)] = -ra * dt;          /* -R[a]x dt */
            Phi[MAT15_IDX(3 + i, 9 + j)] = -R[3 * i + j] * dt;/* -R dt     */
            Phi[MAT15_IDX(6 + i, 6 + j)] += -Sw[3 * i + j] * dt; /* I - [w]x dt */
        }
    }
}

static void cov_update(s4f_model_t *m)
{
    float Phi[ERR_DIM * ERR_DIM];
    float qd[ERR_DIM];
    int i;
    if (m->dt_acc <= 0.0f) return;
    build_phi(m, m->dt_acc, Phi);
    for (i = 0; i < ERR_DIM; i++) qd[i] = m->prm.qc[i / 3] * m->dt_acc;
    m->n_pd_viol += thornton(m->ekf.U, m->ekf.D, Phi, qd);
    m->dt_acc = 0.0f;
    m->cnt = 0;
}

/* Bierman scalar update, H = e_i0 (unit row), see ICD Role 3 */
static void bierman(es_ekf_t *e, int i0, float y, float Rv)
{
    float f[ERR_DIM], g[ERR_DIM], k[ERR_DIM];
    float alpha = Rv, innov;
    int i, j;
    for (j = 0; j < ERR_DIM; j++) {
        f[j] = (j >= i0) ? e->U[U_IDX(i0, j)] : 0.0f;
        g[j] = e->D[j] * f[j];
        k[j] = 0.0f;
    }
    for (j = i0; j < ERR_DIM; j++) {
        float an = alpha + f[j] * g[j];
        float dn = e->D[j] * alpha / an;
        float lam = -f[j] / alpha;
        k[j] = g[j];
        for (i = 0; i < j; i++) {
            float tmp = e->U[U_IDX(i, j)];
            e->U[U_IDX(i, j)] = tmp + lam * k[i];
            k[i] += tmp * g[j];
        }
        e->D[j] = dn;
        alpha = an;
    }
    innov = y - e->dx[i0];
    for (j = 0; j < ERR_DIM; j++) e->dx[j] += k[j] * innov / alpha;
}

static void inject_reset(es_ekf_t *e, float *peak)
{
    float dq[4], qn[4];
    int i;
    for (i = 0; i < ERR_DIM; i++) {
        float a = fabsf(e->dx[i]) / (float)ldexp(1.0, s4_scale_log2[i]);
        if (a > *peak) *peak = a;
    }
    for (i = 0; i < 3; i++) {
        e->x_nom.p[i] += e->dx[i];
        e->x_nom.v[i] += e->dx[3 + i];
        e->x_nom.ba[i] += e->dx[9 + i];
        e->x_nom.bg[i] += e->dx[12 + i];
    }
    dq[0] = 1.0f;
    for (i = 0; i < 3; i++) dq[1 + i] = 0.5f * e->dx[6 + i];
    quat_mul(e->x_nom.q, dq, qn);
    quat_normalize(qn);
    memcpy(e->x_nom.q, qn, sizeof qn);
    for (i = 0; i < ERR_DIM; i++) e->dx[i] = 0.0f;
}

void s4f_step(s4f_model_t *m, const sensor_frame_t *f)
{
    nominal_state_t *x = &m->ekf.x_nom;
    float an[3], dq[4], qn[4];
    float dt = f->dt;
    int i, j;

    /* ---- nominal propagation (ICD Role 1) ---- */
    for (i = 0; i < 3; i++) {
        m->a_unb[i] = f->imu_acc[i] - x->ba[i];
        m->w_unb[i] = f->imu_gyr[i] - x->bg[i];
    }
    quat_to_R(x->q, m->R);
    for (i = 0; i < 3; i++) {
        an[i] = 0.0f;
        for (j = 0; j < 3; j++) an[i] += m->R[3 * i + j] * m->a_unb[j];
    }
    an[2] += (float)S4_GRAVITY_MPS2;
    for (i = 0; i < 3; i++) {
        x->p[i] += x->v[i] * dt + 0.5f * an[i] * dt * dt;
        x->v[i] += an[i] * dt;
    }
    dq[0] = 1.0f;
    for (i = 0; i < 3; i++) dq[1 + i] = 0.5f * m->w_unb[i] * dt;
    quat_mul(x->q, dq, qn);
    quat_normalize(qn);
    memcpy(x->q, qn, sizeof qn);

    /* ---- covariance time update (decimated) ---- */
    m->dt_acc += dt;
    if (++m->cnt >= m->prm.cov_decim) cov_update(m);

    /* ---- measurement update ---- */
    if (f->type == SENS_GNSS_POS_VEL) {
        cov_update(m);
        for (i = 0; i < 3; i++) bierman(&m->ekf, i, f->meas[i] - x->p[i], m->prm.r_pos);
        for (i = 0; i < 3; i++) bierman(&m->ekf, 3 + i, f->meas[3 + i] - x->v[i], m->prm.r_vel);
        inject_reset(&m->ekf, &m->peak_dx_n);
    } else if (f->type == SENS_BARO) {
        cov_update(m);
        bierman(&m->ekf, 2, f->meas[0] - x->p[2], m->prm.r_baro);
        inject_reset(&m->ekf, &m->peak_dx_n);
    }
    /* SENS_MAG: PENDING (not implemented in System 4 study pipeline) */
}
