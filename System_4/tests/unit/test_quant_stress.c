/* Test 5 -- quantization stress: tiny, between-LSB, near zero, max, min, per format. */
#include <stdint.h>
#include <stdio.h>
#include "ekf_fixed_math.h"
#include "test_util.h"

typedef struct { const char *n; uint8_t wl, frac; } fmt_t;

int main(void)
{
    static const fmt_t F[] = {{"Q1.15", 16, 15}, {"Q2.14", 16, 14}, {"Q1.23", 24, 23},
                              {"Q2.22", 24, 22}, {"Q2.30", 32, 30}, {"Q16.16", 32, 16}};
    size_t f;
    int i;
    TEST_BEGIN("quantization stress");
    for (f = 0; f < sizeof F / sizeof F[0]; f++) {
        const fmt_t *q = &F[f];
        double lsb = ldexp(1.0, -(int)q->frac);
        double maxv = (double)fxp_wl_max(q->wl) * lsb;
        double minv = (double)fxp_wl_min(q->wl) * lsb;
        double maxerr = 0;
        printf("  %s: lsb %.3e range [%g, %g]\n", q->n, lsb, minv, maxv);

        /* tiny values vanish, but never produce a spurious non-zero */
        CHECK_EQ(fxp_from_double_q(lsb * 1e-6, q->wl, q->frac), 0, "far below 1/2 LSB -> 0");
        CHECK_EQ(fxp_from_double_q(-lsb * 1e-6, q->wl, q->frac), 0, "negative tiny -> 0 (no -0 bias)");
        CHECK_EQ(fxp_from_double_q(lsb * 0.49, q->wl, q->frac), 0, "0.49 LSB -> 0");
        CHECK_EQ(fxp_from_double_q(lsb * 0.51, q->wl, q->frac), 1, "0.51 LSB -> 1");
        CHECK_EQ(fxp_from_double_q(-lsb * 0.51, q->wl, q->frac), -1, "-0.51 LSB -> -1");
        CHECK_EQ(fxp_from_double_q(0.0, q->wl, q->frac), 0, "zero");
        /* exactly between two LSBs: half away from zero (symmetric) */
        CHECK_EQ(fxp_from_double_q(lsb * 1.5, q->wl, q->frac), 2, "1.5 LSB -> 2");
        CHECK_EQ(fxp_from_double_q(-lsb * 1.5, q->wl, q->frac), -2, "-1.5 LSB -> -2");
        CHECK_EQ(fxp_from_double_q(lsb * 2.5, q->wl, q->frac), 3, "2.5 LSB -> 3 (half away, not banker's)");
        /* extremes */
        CHECK_EQ((long long)fxp_from_double_q(maxv, q->wl, q->frac), fxp_wl_max(q->wl), "max exact");
        CHECK_EQ((long long)fxp_from_double_q(minv, q->wl, q->frac), fxp_wl_min(q->wl), "min exact");
        CHECK_EQ((long long)fxp_from_double_q(maxv * 2.0 + 1.0, q->wl, q->frac), fxp_wl_max(q->wl), "beyond max saturates");
        CHECK_EQ((long long)fxp_from_double_q(minv * 2.0 - 1.0, q->wl, q->frac), fxp_wl_min(q->wl), "beyond min saturates");
        CHECK_EQ((long long)fxp_from_double_q(1e300, q->wl, q->frac), fxp_wl_max(q->wl), "huge saturates (no UB)");
        CHECK_EQ((long long)fxp_from_double_q(-1e300, q->wl, q->frac), fxp_wl_min(q->wl), "-huge saturates (no UB)");
        /* error bound over random in-range values */
        for (i = 0; i < 100000; i++) {
            double x = minv + (maxv - minv) * ((double)(t_rand() >> 8) / 16777216.0);
            double e = fabs(fxp_to_double_q(fxp_from_double_q(x, q->wl, q->frac), q->frac) - x);
            if (e > maxerr) maxerr = e;
        }
        CHECK(maxerr <= 0.5 * lsb * (1.0 + 1e-9), "max rounding error <= 1/2 LSB");
        CHECK(maxerr >= 0.45 * lsb, "bound is nearly attained (histogram uniform)");
        printf("    measured max |error| = %.4f LSB\n", maxerr / lsb);
    }
    /* saturation of a normalized +1.0 differs per format: documented behaviour */
    CHECK_EQ(fxp_from_double_q(1.0, 16, 15), 32767, "Q1.15: +1.0 saturates to 1-LSB");
    CHECK_EQ(fxp_from_double_q(1.0, 24, 23), 8388607, "Q1.23: +1.0 saturates to 1-LSB");
    CHECK_EQ(fxp_from_double_q(1.0, 32, 30), 1 << 30, "Q2.30: +1.0 exact");
    return TEST_END();
}
