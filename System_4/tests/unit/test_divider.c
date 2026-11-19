/* Divider tests: baseline vs restoring vs non-restoring, errors, div0, saturation, latency. */
#include <stdint.h>
#include <stdio.h>
#include "ekf_fixed_divider.h"
#include "ekf_fixed_math.h"
#include "test_util.h"

int main(void)
{
    const div_mode_t modes[3] = {DIV_BASELINE, DIV_RESTORING, DIV_NONRESTORING};
    int m, i;
    TEST_BEGIN("ekf_fixed_divider");

    /* unsigned cores equal native division (including remainder) */
    for (i = 0; i < 20000; i++) {
        uint64_t n = ((uint64_t)t_rand() << 32) | t_rand();
        uint64_t d = (((uint64_t)t_rand() << 32) | t_rand()) >> (1 + (t_rand() % 40)); /* < 2^63 */
        uint64_t r1, r2;
        if (d == 0) d = 1;
        if (d >= ((uint64_t)1 << 62)) d >>= 2; /* non-restoring keeps |r| < 2d < 2^63 */
        {
            uint64_t q1 = divu_restoring(n, d, 64, &r1);
            uint64_t q2 = divu_nonrestoring(n, d, 64, &r2);
            if (q1 != n / d || r1 != n % d) { CHECK(0, "restoring != native"); break; }
            if (q2 != n / d || r2 != n % d) { CHECK(0, "non-restoring != native"); break; }
        }
    }
    CHECK(1, "20k random 64-bit divu checks");

    CHECK_EQ(div_bit_width(0), 0, "bit width 0");
    CHECK_EQ(div_bit_width(1), 1, "bit width 1");
    CHECK_EQ(div_bit_width(255), 8, "bit width 255");

    for (m = 0; m < 3; m++) {
        div_info_t inf;
        int mism = 0;
        double maxe = 0;
        printf("  mode %s\n", div_mode_name(modes[m]));
        /* known values (Q16.16) */
        CHECK_EQ(div_q(3 << 16, 2 << 16, 16, 16, 16, modes[m], &inf), 3 << 15, "3/2");
        CHECK_EQ(div_q(-3 * 65536, 2 << 16, 16, 16, 16, modes[m], &inf), -(3 << 15), "-3/2");
        CHECK_EQ(div_q(3 << 16, -2 * 65536, 16, 16, 16, modes[m], &inf), -(3 << 15), "3/-2");
        CHECK_EQ(div_q(-3 * 65536, -2 * 65536, 16, 16, 16, modes[m], &inf), 3 << 15, "-3/-2");
        CHECK_EQ(div_q(0, 7 << 16, 16, 16, 16, modes[m], &inf), 0, "0/x");
        CHECK_EQ(div_q(1 << 16, 3 << 16, 16, 16, 16, modes[m], &inf), 21845, "1/3");
        CHECK_EQ(div_q(-1 * 65536, 3 << 16, 16, 16, 16, modes[m], &inf), -21845, "-1/3 toward zero");
        /* edge: INT32_MIN / -1 saturates */
        CHECK_EQ(div_q(INT32_MIN, -1, 0, 0, 0, modes[m], &inf), INT32_MAX, "INT32_MIN/-1 saturates");
        CHECK_EQ(div_q(INT32_MAX, 1, 16, 16, 16, modes[m], &inf), INT32_MAX, "max/1 (Q16.16: 2^31/..) ok or sat");
        CHECK_EQ(div_q(INT32_MAX, 1 << 14, 16, 16, 16, modes[m], &inf), INT32_MAX, "max/0.25 saturates");
        /* division by zero */
        fxp_stats_reset();
        CHECK_EQ(div_q(5 << 16, 0, 16, 16, 16, modes[m], &inf), INT32_MAX, "x/0 -> +max");
        CHECK_EQ(inf.div0, 1, "div0 flag set");
        CHECK_EQ(div_q(-5 * 65536, 0, 16, 16, 16, modes[m], &inf), INT32_MIN, "-x/0 -> min");
        CHECK_EQ(div_q(0, 0, 16, 16, 16, modes[m], &inf), 0, "0/0 -> 0");
        CHECK_EQ(fxp_stats_total().div0, 3, "stats counted div0");
        /* random: bit-exact against baseline fxp_div_q, <=1 LSB from exact */
        for (i = 0; i < 50000; i++) {
            int32_t a = (int32_t)(t_rand() >> 8) - (1 << 23);
            int32_t b = (int32_t)(t_rand() >> 10) - (1 << 21);
            int32_t r, rb;
            if (b == 0) b = 1;
            r = div_q(a, b, 16, 16, 16, modes[m], &inf);
            rb = fxp_div_q(a, b, 16, 16, 16);
            if (r != rb) mism++;
            if (r != INT32_MAX && r != INT32_MIN) {
                double e = fabs((double)r - ((double)a / (double)b) * 65536.0);
                if (e > maxe) maxe = e;
            }
        }
        CHECK_EQ(mism, 0, "bit-exact with baseline");
        CHECK(maxe <= 1.0 + 1e-9, "error <= 1 LSB (truncation)");
        printf("    max error %.4f LSB; one sample: iterations=%u cycles_est=%u\n", maxe,
               inf.iterations, inf.cycles_est);
    }

    /* rounding variant */
    CHECK_EQ((long long)divs64_round(7, 2, DIV_RESTORING, (div_info_t *)0), 4, "round(3.5)=4");
    CHECK_EQ((long long)divs64_round(-7, 2, DIV_NONRESTORING, (div_info_t *)0), -4, "round(-3.5)=-4");
    CHECK_EQ((long long)divs64_round(5, 3, DIV_BASELINE, (div_info_t *)0), 2, "round(1.67)=2");
    CHECK_EQ((long long)divs64_round(4, 3, DIV_RESTORING, (div_info_t *)0), 1, "round(1.33)=1");

    /* iteration counts follow dividend width (fixed-latency radix-2) */
    {
        div_info_t a, b;
        (void)divs64(1000, 7, DIV_RESTORING, &a);
        (void)divs64(1000000000LL, 7, DIV_RESTORING, &b);
        CHECK(b.iterations > a.iterations, "wider dividend -> more iterations");
        CHECK_EQ(a.iterations, div_bit_width(1000), "iterations == dividend bit width");
        (void)divs64(1000, 7, DIV_NONRESTORING, &b);
        CHECK_EQ(b.cycles_est, a.cycles_est + 1, "non-restoring has +1 correction cycle");
    }
    return TEST_END();
}
