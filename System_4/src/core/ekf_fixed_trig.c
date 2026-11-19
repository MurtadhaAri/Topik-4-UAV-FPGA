/*
 * ekf_fixed_trig.c -- CORDIC sin/cos.  Core uses only int32 add/sub/shift and
 * a lookup table; no sin()/cos()/atan() library calls.
 */
#include "ekf_fixed_trig.h"

#include "ekf_fixed_math.h"

/* atan(2^-i) in Q3.29 (round to nearest), i = 0..31.  Entries >= 29 are
 * <= 1 LSB, which is why CORDIC_MAX_ITERATIONS is 30. */
const int32_t CORDIC_ATAN_Q29[32] = {
    421657428, 248918915, 131521918, 66762579, 33510843, 16771758, 8387925, 4194219,
    2097141,   1048575,   524288,    262144,   131072,   65536,    32768,   16384,
    8192,      4096,      2048,      1024,     512,      256,      128,     64,
    32,        16,        8,         4,        2,        1,        1,       0};

void cordic_sincos_n(int32_t angle, uint8_t iterations, int32_t *sin_out, int32_t *cos_out)
{
    int64_t a = (int64_t)angle;
    int32_t x, y, z;
    int neg_cos = 0;
    uint8_t i;

    if (iterations > CORDIC_MAX_ITERATIONS) iterations = (uint8_t)CORDIC_MAX_ITERATIONS;

    /* wrap into [-pi, pi]; |angle| < 4 < 2*pi hence at most one wrap */
    if (a > (int64_t)CORDIC_PI_Q29) a -= CORDIC_TWOPI_Q29;
    else if (a < -(int64_t)CORDIC_PI_Q29) a += CORDIC_TWOPI_Q29;

    /* quadrant reduction to [-pi/2, pi/2] */
    if (a > (int64_t)CORDIC_HALFPI_Q29) {
        a = (int64_t)CORDIC_PI_Q29 - a;
        neg_cos = 1;
    } else if (a < -(int64_t)CORDIC_HALFPI_Q29) {
        a = -(int64_t)CORDIC_PI_Q29 - a;
        neg_cos = 1;
    }

    x = CORDIC_K_Q30;
    y = 0;
    z = (int32_t)a;

    /* fixed iteration count: the loop bound is the (constant-per-call) N */
    for (i = 0; i < iterations; i++) {
        int32_t xs = x >> i; /* arithmetic shift of int32 (GCC, documented) */
        int32_t ys = y >> i;
        if (z >= 0) {
            x -= ys;
            y += xs;
            z -= CORDIC_ATAN_Q29[i];
        } else {
            x += ys;
            y -= xs;
            z += CORDIC_ATAN_Q29[i];
        }
    }

    if (neg_cos) x = -x; /* |x| <= 2^30 : cannot overflow */
    if (sin_out) *sin_out = y;
    if (cos_out) *cos_out = x;
}

void cordic_sincos(int32_t angle, int32_t *sin_out, int32_t *cos_out)
{
    cordic_sincos_n(angle, (uint8_t)CORDIC_ITERATIONS, sin_out, cos_out);
}

/* BEGIN-CONVERSION */
int32_t cordic_angle_from_double(double rad)
{
    return fxp_from_double_q(rad, 32u, (uint8_t)CORDIC_ANGLE_FRAC);
}
/* END-CONVERSION */
