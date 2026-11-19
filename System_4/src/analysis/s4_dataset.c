/*
 * s4_dataset.c -- synthetic truth + sensor generator (double precision; this
 * is analysis tooling, not part of the fixed-point core).
 */
#include "s4_dataset.h"

#include <math.h>
#include <string.h>

#include "ekf_fixed_config.h"

/* ---- deterministic PRNG (xorshift32) + Box-Muller ---- */
typedef struct { uint32_t s; int have; double spare; } rng_t;

static double rng_u(rng_t *r)
{
    r->s ^= r->s << 13;
    r->s ^= r->s >> 17;
    r->s ^= r->s << 5;
    return ((double)r->s + 0.5) / 4294967296.0;
}

static double rng_n(rng_t *r)
{
    double u1, u2, m;
    if (r->have) {
        r->have = 0;
        return r->spare;
    }
    u1 = rng_u(r);
    u2 = rng_u(r);
    m = sqrt(-2.0 * log(u1));
    r->spare = m * sin(6.283185307179586 * u2);
    r->have = 1;
    return m * cos(6.283185307179586 * u2);
}

const char *s4_scenario_name(s4_scenario_t sc)
{
    switch (sc) {
    case S4_SC_STATIC: return "static";
    case S4_SC_CONST_VEL: return "constvel";
    case S4_SC_DYNAMIC: return "dynamic";
    case S4_SC_HIGH_RANGE: return "highrange";
    default: return "?";
    }
}

static void qmul(const double *p, const double *q, double *r)
{
    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
    r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
    r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
    r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

static void qexp(const double *phi, double *dq)
{
    double a = sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
    double s = (a > 1e-12) ? sin(0.5 * a) / a : 0.5;
    dq[0] = cos(0.5 * a);
    dq[1] = s * phi[0];
    dq[2] = s * phi[1];
    dq[3] = s * phi[2];
}

static void q_to_R(const double *q, double *R)
{
    double w = q[0], x = q[1], y = q[2], z = q[3];
    R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - w * z); R[2] = 2 * (x * z + w * y);
    R[3] = 2 * (x * y + w * z); R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - w * x);
    R[6] = 2 * (x * z - w * y); R[7] = 2 * (y * z + w * x); R[8] = 1 - 2 * (x * x + y * y);
}

void s4_dataset_filter_init(s4_scenario_t sc, float v0[3])
{
    v0[0] = v0[1] = v0[2] = 0.0f;
    if (sc == S4_SC_CONST_VEL) {
        v0[0] = 10.0f; v0[1] = 5.0f; v0[2] = -0.5f;
    } else if (sc == S4_SC_DYNAMIC) {
        v0[0] = 5.0f;
    }
}

int s4_dataset_generate(s4_scenario_t sc, double duration_s, sensor_frame_t *frames,
                        s4_truth_t *truth, int max_frames)
{
    const double dt = S4_IMU_DT;
    const double g = S4_GRAVITY_MPS2;
    int n = (int)(duration_s / dt + 0.5);
    int noise = (sc != S4_SC_STATIC);
    double p[3] = {0, 0, 0}, v[3] = {0, 0, 0}, q[4] = {1, 0, 0, 0};
    double ba[3] = {0, 0, 0}, bg[3] = {0, 0, 0};
    double acc_amp[3] = {0, 0, 0}, w_amp[3] = {0, 0, 0}, w_off[3] = {0, 0, 0};
    rng_t rng;
    int k, i;

    rng.s = 0x1234ABCDu + (uint32_t)sc * 7919u;
    rng.have = 0;
    rng.spare = 0.0;
    if (n > max_frames) n = max_frames;

    if (sc == S4_SC_CONST_VEL) {
        v[0] = 10.0; v[1] = 5.0; v[2] = -0.5;
    } else if (sc == S4_SC_DYNAMIC || sc == S4_SC_HIGH_RANGE) {
        v[0] = 5.0;
        acc_amp[0] = 1.5; acc_amp[1] = 1.0; acc_amp[2] = 0.5;
        w_amp[0] = 0.4; w_amp[1] = 0.3; w_amp[2] = 0.2; w_off[2] = 0.1;
        ba[0] = 0.05; ba[1] = -0.03; ba[2] = 0.02;
        bg[0] = 0.002; bg[1] = -0.001; bg[2] = 0.003;
    }
    if (sc == S4_SC_HIGH_RANGE) { /* ~95% of the ICD error bounds, filter starts at 0 */
        double rv[3] = {0.28, 0.28, 0.28}, dq[4], t[4];
        p[0] = 9.5; p[1] = -9.5; p[2] = 9.5;
        v[0] += 9.5; v[1] = -9.5; v[2] = 9.5;
        qexp(rv, dq);
        qmul(q, dq, t);
        memcpy(q, t, sizeof t);
        ba[0] = 1.9; ba[1] = -1.9; ba[2] = 1.9;
        bg[0] = 0.19; bg[1] = -0.19; bg[2] = 0.19;
    }

    for (k = 1; k <= n; k++) {
        double tm = (k - 0.5) * dt;
        double an[3], w[3], phi[3], dq[4], qn[4], qm[4], R[9], fn[3], fb[3];
        sensor_frame_t *f = &frames[k - 1];

        an[0] = acc_amp[0] * sin(0.4 * tm);
        an[1] = acc_amp[1] * cos(0.3 * tm);
        an[2] = acc_amp[2] * sin(0.2 * tm);
        w[0] = w_amp[0] * sin(0.7 * tm) + w_off[0];
        w[1] = w_amp[1] * cos(0.5 * tm) + w_off[1];
        w[2] = w_amp[2] * sin(0.3 * tm) + w_off[2];

        for (i = 0; i < 3; i++) phi[i] = 0.5 * w[i] * dt;
        qexp(phi, dq);
        qmul(q, dq, qm); /* attitude at mid-interval */
        for (i = 0; i < 3; i++) phi[i] = w[i] * dt;
        qexp(phi, dq);
        qmul(q, dq, qn);

        q_to_R(qm, R);
        for (i = 0; i < 3; i++) fn[i] = an[i] - ((i == 2) ? g : 0.0); /* a_n - g_n */
        for (i = 0; i < 3; i++) fb[i] = R[i] * fn[0] + R[3 + i] * fn[1] + R[6 + i] * fn[2]; /* R^T */

        memset(f, 0, sizeof *f);
        f->timestamp_us = (uint32_t)k * 2000u;
        f->dt = (float)dt;
        f->type = SENS_IMU_ONLY;
        for (i = 0; i < 3; i++) {
            f->imu_acc[i] = (float)(fb[i] + ba[i] + (noise ? 0.05 * rng_n(&rng) : 0.0));
            f->imu_gyr[i] = (float)(w[i] + bg[i] + (noise ? 0.003 * rng_n(&rng) : 0.0));
            p[i] += v[i] * dt + 0.5 * an[i] * dt * dt;
            v[i] += an[i] * dt;
        }
        memcpy(q, qn, sizeof qn);

        if (k % 100 == 0) {
            f->type = SENS_GNSS_POS_VEL;
            for (i = 0; i < 3; i++) {
                f->meas[i] = (float)(p[i] + (noise ? 0.5 * rng_n(&rng) : 0.0));
                f->meas[3 + i] = (float)(v[i] + (noise ? 0.2 * rng_n(&rng) : 0.0));
            }
        } else if (k % 25 == 0) {
            f->type = SENS_BARO;
            f->meas[0] = (float)(p[2] + (noise ? 0.4 * rng_n(&rng) : 0.0)); /* z = -altitude = pD */
        }

        for (i = 0; i < 3; i++) {
            truth[k - 1].p[i] = p[i];
            truth[k - 1].v[i] = v[i];
        }
        for (i = 0; i < 4; i++) truth[k - 1].q[i] = q[i];
    }
    return n;
}
