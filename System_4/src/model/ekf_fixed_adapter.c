/*
 * ekf_fixed_adapter.c -- float <-> fixed conversion with normalization.
 * Uses double: allowed here (I/O + analysis only).
 */
#include "ekf_fixed_adapter.h"

#include <math.h>

#include "ekf_fixed_math.h"

const int8_t s4_scale_log2[ERR_DIM] = {4, 4, 4, 4, 4, 4, -1, -1, -1, 1, 1, 1, -2, -2, -2};

double s4_scale(int i) { return ldexp(1.0, s4_scale_log2[i]); }

/* ------------------------------------------------ configuration table */
/* PROPOSED DESIGN PARAMETERS (nominal path identical in all three: ICD Role 4
 * suggests Q16.16 position and Q2.30 quaternion).  Integer-bit choices of the
 * error path come from results/csv/range_analysis.csv (measured, 60 s runs):
 *   dx  Q4.x : peak |dx_n| 3.5 (highrange)         -> 3 bits needed, 4 used
 *   U   Q7.x : max |U_n| 30.0 (yaw unobserved, no magnetometer) -> 6 needed, 7 used
 *   D   Q1.x : max D_n 0.21                        -> 1 bit needed
 *   Phi Q1.x : |dPhi_n| << 1                                                */
static const ekf_qcfg_t k_cfg[S4_NUM_CFG] = {
    /* name   p        v        q        b        dx       U         D        Phi    */
    {"W16", {32, 16}, {32, 24}, {32, 30}, {32, 28}, {16, 12}, {16, 9}, {16, 15}, {16, 15}},
    {"W24", {32, 16}, {32, 24}, {32, 30}, {32, 28}, {24, 20}, {24, 17}, {24, 23}, {24, 23}},
    {"W32", {32, 16}, {32, 24}, {32, 30}, {32, 28}, {32, 28}, {32, 25}, {32, 31}, {32, 31}},
};

const ekf_qcfg_t *s4_cfg_get(int index)
{
    return (index >= 0 && index < S4_NUM_CFG) ? &k_cfg[index] : (const ekf_qcfg_t *)0;
}

/* --------------------------------------------------------- nominal state */
static void nom_to_fixed(const nominal_state_t *s, nominal_state_fixed_t *d, const ekf_qcfg_t *c)
{
    int i;
    for (i = 0; i < 3; i++) {
        d->p[i] = fxp_from_double_q((double)s->p[i], c->p.wl, c->p.frac);
        d->v[i] = fxp_from_double_q((double)s->v[i], c->v.wl, c->v.frac);
        d->ba[i] = fxp_from_double_q((double)s->ba[i], c->b.wl, c->b.frac);
        d->bg[i] = fxp_from_double_q((double)s->bg[i], c->b.wl, c->b.frac);
    }
    for (i = 0; i < 4; i++) d->q[i] = fxp_from_double_q((double)s->q[i], c->q.wl, c->q.frac);
}

static void nom_to_float(const nominal_state_fixed_t *s, nominal_state_t *d, const ekf_qcfg_t *c)
{
    int i;
    for (i = 0; i < 3; i++) {
        d->p[i] = (float)fxp_to_double_q(s->p[i], c->p.frac);
        d->v[i] = (float)fxp_to_double_q(s->v[i], c->v.frac);
        d->ba[i] = (float)fxp_to_double_q(s->ba[i], c->b.frac);
        d->bg[i] = (float)fxp_to_double_q(s->bg[i], c->b.frac);
    }
    for (i = 0; i < 4; i++) d->q[i] = (float)fxp_to_double_q(s->q[i], c->q.frac);
}

/* ------------------------------------------------------------- full state */
void ekf_float_to_fixed(const es_ekf_t *src, es_ekf_fixed_t *dst, const ekf_qcfg_t *c)
{
    int i, j;
    nom_to_fixed(&src->x_nom, &dst->x_nom, c);
    for (i = 0; i < ERR_DIM; i++) {
        dst->dx[i] = fxp_from_double_q(ldexp((double)src->dx[i], -s4_scale_log2[i]), c->dx.wl,
                                       c->dx.frac);
        dst->D[i] = fxp_from_double_q(ldexp((double)src->D[i], -2 * s4_scale_log2[i]), c->D.wl,
                                      c->D.frac);
    }
    for (j = 0; j < ERR_DIM; j++)
        for (i = 0; i <= j; i++)
            dst->U[U_IDX(i, j)] = fxp_from_double_q(
                ldexp((double)src->U[U_IDX(i, j)], s4_scale_log2[j] - s4_scale_log2[i]), c->U.wl,
                c->U.frac);
}

void ekf_fixed_to_float(const es_ekf_fixed_t *src, es_ekf_t *dst, const ekf_qcfg_t *c)
{
    int i, j;
    nom_to_float(&src->x_nom, &dst->x_nom, c);
    for (i = 0; i < ERR_DIM; i++) {
        dst->dx[i] = (float)ldexp(fxp_to_double_q(src->dx[i], c->dx.frac), s4_scale_log2[i]);
        dst->D[i] = (float)ldexp(fxp_to_double_q(src->D[i], c->D.frac), 2 * s4_scale_log2[i]);
    }
    for (j = 0; j < ERR_DIM; j++)
        for (i = 0; i <= j; i++)
            dst->U[U_IDX(i, j)] = (float)ldexp(fxp_to_double_q(src->U[U_IDX(i, j)], c->U.frac),
                                               s4_scale_log2[i] - s4_scale_log2[j]);
}

/* ---------------------------------------------------------- sensor packet */
void sensor_frame_to_fixed(const sensor_frame_t *s, sensor_frame_fixed_t *d, const ekf_qcfg_t *c)
{
    int i;
    d->timestamp_us = s->timestamp_us;
    d->type = s->type;
    d->dt = fxp_from_double_q((double)s->dt, 32u, (uint8_t)S4_DT_FRAC);
    for (i = 0; i < 3; i++) {
        d->imu_acc[i] = fxp_from_double_q((double)s->imu_acc[i], 32u, (uint8_t)S4_IMU_FRAC);
        d->imu_gyr[i] = fxp_from_double_q((double)s->imu_gyr[i], 32u, (uint8_t)S4_IMU_FRAC);
        d->meas[i] = fxp_from_double_q((double)s->meas[i], c->p.wl, c->p.frac);
        d->meas[3 + i] = fxp_from_double_q((double)s->meas[3 + i], c->v.wl, c->v.frac);
    }
}

/* ------------------------------------------------ normalized access helpers */
double s4_float_Un(const es_ekf_t *e, int i, int j)
{
    return ldexp((double)e->U[U_IDX(i, j)], s4_scale_log2[j] - s4_scale_log2[i]);
}

double s4_float_Dn(const es_ekf_t *e, int i)
{
    return ldexp((double)e->D[i], -2 * s4_scale_log2[i]);
}

double s4_fixed_Un(const es_ekf_fixed_t *e, const ekf_qcfg_t *c, int i, int j)
{
    return fxp_to_double_q(e->U[U_IDX(i, j)], c->U.frac);
}

double s4_fixed_Dn(const es_ekf_fixed_t *e, const ekf_qcfg_t *c, int i)
{
    return fxp_to_double_q(e->D[i], c->D.frac);
}
