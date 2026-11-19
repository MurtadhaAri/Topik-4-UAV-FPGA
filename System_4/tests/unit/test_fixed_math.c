/* Test 6 -- arithmetic unit tests of ekf_fixed_math (default Q16.16 unless stated). */
#include <stdint.h>
#include "ekf_fixed_math.h"
#include "test_util.h"

#define ONE (1 << 16)

int main(void)
{
    int i;
    fxp_counters_t c;
    TEST_BEGIN("ekf_fixed_math");

    /* ---------------- saturation primitive ---------------- */
    fxp_stats_reset();
    CHECK_EQ(fxp_saturate(0), 0, "sat 0");
    CHECK_EQ(fxp_saturate(INT32_MAX), INT32_MAX, "sat max");
    CHECK_EQ(fxp_saturate((int64_t)INT32_MAX + 1), INT32_MAX, "sat max+1");
    CHECK_EQ(fxp_saturate((int64_t)INT32_MIN - 1), INT32_MIN, "sat min-1");
    CHECK_EQ(fxp_saturate(INT64_MAX), INT32_MAX, "sat int64 max");
    CHECK_EQ(fxp_saturate(INT64_MIN), INT32_MIN, "sat int64 min");
    c = fxp_stats_total();
    CHECK_EQ(c.sat, 4, "4 saturations counted");
    CHECK_EQ(c.ovf, 4, "4 container overflows counted");
    CHECK_EQ(fxp_sat_wl(40000, 16), 32767, "wl16 clip high");
    CHECK_EQ(fxp_sat_wl(-40000, 16), -32768, "wl16 clip low");
    CHECK_EQ(fxp_sat_wl(32767, 16), 32767, "wl16 edge");
    CHECK_EQ(fxp_sat_wl(8388608, 24), 8388607, "wl24 clip high");
    CHECK_EQ(fxp_sat_wl(-8388609, 24), -8388608, "wl24 clip low");
    fxp_stats_reset();
    (void)fxp_sat_wl(40000, 16);
    c = fxp_stats_total();
    CHECK_EQ(c.sat, 1, "wl16 clip counted as saturation");
    CHECK_EQ(c.ovf, 0, "wl16 clip is NOT a container overflow");

    /* ---------------- add / sub ---------------- */
    CHECK_EQ(fxp_add_sat(5, 7), 12, "add pos");
    CHECK_EQ(fxp_add_sat(-5, -7), -12, "add neg");
    CHECK_EQ(fxp_add_sat(0, 0), 0, "add zero");
    CHECK_EQ(fxp_add_sat(5, -5), 0, "add cancel");
    CHECK_EQ(fxp_add_sat(INT32_MAX - 1, 1), INT32_MAX, "add near max exact");
    CHECK_EQ(fxp_add_sat(INT32_MAX, 1), INT32_MAX, "add overflow saturates");
    CHECK_EQ(fxp_add_sat(INT32_MIN, -1), INT32_MIN, "add underflow saturates");
    CHECK_EQ(fxp_add_sat(INT32_MIN + 1, -1), INT32_MIN, "add near min exact");
    CHECK_EQ(fxp_sub_sat(5, 7), -2, "sub");
    CHECK_EQ(fxp_sub_sat(0, INT32_MIN), INT32_MAX, "sub 0-min saturates");
    CHECK_EQ(fxp_sub_sat(INT32_MIN, 1), INT32_MIN, "sub min-1 saturates");
    CHECK_EQ(fxp_sub_sat(INT32_MAX, -1), INT32_MAX, "sub max+1 saturates");
    CHECK_EQ(fxp_add_sat_wl(30000, 30000, 16), 32767, "add wl16 saturates");

    /* ---------------- mul (floor) and mul_round ---------------- */
    CHECK_EQ(fxp_mul(3 * ONE / 2, 2 * ONE), 3 * ONE, "1.5*2=3");
    CHECK_EQ(fxp_mul(-3 * ONE / 2, 2 * ONE), -3 * ONE, "-1.5*2=-3");
    CHECK_EQ(fxp_mul(0, 12345), 0, "0*x");
    CHECK_EQ(fxp_mul(ONE / 2, ONE / 2), ONE / 4, "0.5*0.5");
    CHECK_EQ(fxp_mul(ONE, 777), 777, "1.0*x identity");
    /* sub-LSB behaviour: documented floor vs half-away rounding */
    CHECK_EQ(fxp_mul(ONE / 2, 1), 0, "floor: +0.5 LSB -> 0");
    CHECK_EQ(fxp_mul_round(ONE / 2, 1), 1, "round: +0.5 LSB -> +1 (half away)");
    CHECK_EQ(fxp_mul(-ONE / 2, 1), -1, "floor: -0.5 LSB -> -1 (bias to -inf)");
    CHECK_EQ(fxp_mul_round(-ONE / 2, 1), -1, "round: -0.5 LSB -> -1 (symmetric)");
    CHECK_EQ(fxp_mul_round(ONE / 4, 1), 0, "round: +0.25 LSB -> 0");
    CHECK_EQ(fxp_mul_round(-ONE / 4, 1), 0, "round: -0.25 LSB -> 0");
    /* large / saturating */
    CHECK_EQ(fxp_mul(INT32_MAX, ONE), INT32_MAX, "max*1.0 exact");
    CHECK_EQ(fxp_mul(INT32_MAX, 2 * ONE), INT32_MAX, "max*2 saturates");
    CHECK_EQ(fxp_mul(INT32_MIN, 2 * ONE), INT32_MIN, "min*2 saturates");
    CHECK_EQ(fxp_mul(INT32_MIN, -ONE), INT32_MAX, "min*(-1) saturates high");
    CHECK_EQ(fxp_mul(INT32_MIN, ONE), INT32_MIN, "min*1.0 exact");
    CHECK_EQ(fxp_mul_round(INT32_MIN, -ONE), INT32_MAX, "round min*(-1) saturates high");
    CHECK_EQ(MUL_FXP(2 * ONE, 3 * ONE), 6 * ONE, "MUL_FXP macro");
    /* underflow counter: non-zero operands, exact zero result */
    fxp_stats_reset();
    (void)fxp_mul(3, 3);
    c = fxp_stats_total();
    CHECK_EQ(c.underflow, 1, "underflow counted for 3*3 raw in Q16.16");
    /* mixed formats: Q8.24 * Q2.30 -> Q16.16 */
    CHECK_EQ(fxp_mul_q(3 << 23, 1 << 29, 24, 30, 16), (3 << 14), "1.5 (Q8.24) * 0.5 (Q2.30) -> 0.75 Q16.16");
    CHECK_EQ(fxp_mul_round_q(1 << 24, 1 << 30, 24, 30, 30), 1 << 30, "Q8.24 1.0 * Q2.30 1.0 -> Q2.30 1.0");
    CHECK_EQ(fxp_mul_round_wl(20000, 20000, 15, 15, 15, 16), 12207, "wl16 Q1.15 mul rounds (0.61*0.61=0.372)");
    CHECK_EQ(fxp_mul_round_wl(32767, 32767, 0, 0, 0, 16), 32767, "wl16 product saturates");

    /* ---------------- division ---------------- */
    CHECK_EQ(fxp_div(3 * ONE, 2 * ONE), 3 * ONE / 2, "3/2");
    CHECK_EQ(fxp_div(-3 * ONE, 2 * ONE), -3 * ONE / 2, "-3/2");
    CHECK_EQ(fxp_div(0, 5 * ONE), 0, "0/x");
    CHECK_EQ(fxp_div(ONE, 3 * ONE), 21845, "1/3 truncated");
    CHECK_EQ(fxp_div(-ONE, 3 * ONE), -21845, "-1/3 truncated toward zero");
    CHECK_EQ(fxp_div(-1, 2 * ONE), 0, "tiny negative / 2.0 -> 0 (toward zero)");
    fxp_stats_reset();
    CHECK_EQ(fxp_div(ONE, 0), INT32_MAX, "x/0 -> +max");
    CHECK_EQ(fxp_div(-ONE, 0), INT32_MIN, "-x/0 -> min");
    CHECK_EQ(fxp_div(0, 0), 0, "0/0 -> 0");
    c = fxp_stats_total();
    CHECK_EQ(c.div0, 3, "three div-by-zero events counted");
    CHECK_EQ(fxp_div(INT32_MAX, ONE / 2), INT32_MAX, "max/0.5 saturates");
    CHECK_EQ(fxp_div_q(1 << 24, 3 << 24, 24, 24, 30), 357913941, "Q8.24 1/3 -> Q2.30");
    CHECK_EQ(fxp_div_q(5, 1, 40, 0, 0), 0, "quotient below 2^-32 of target LSB -> 0");
    CHECK_EQ(fxp_div_q(5, 1, 0, 31, 0), INT32_MAX, "5<<31 / 1 saturates");

    /* ---------------- arithmetic shift ---------------- */
    CHECK_EQ(fxp_shift_right(8, 1), 4, "8>>1");
    CHECK_EQ(fxp_shift_right(-8, 1), -4, "-8>>1");
    CHECK_EQ(fxp_shift_right(-5, 1), -3, "floor(-2.5)=-3");
    CHECK_EQ(fxp_shift_right(-1, 31), -1, "-1>>31");
    CHECK_EQ(fxp_shift_right(-1, 40), -1, "-1>>40 -> -1");
    CHECK_EQ(fxp_shift_right(12345, 40), 0, "x>>40 -> 0");
    CHECK_EQ(fxp_shift_right(INT32_MIN, 31), -1, "min>>31");
    CHECK_EQ((long long)fxp_rshift_round(5, 1), 3, "round(2.5)=3");
    CHECK_EQ((long long)fxp_rshift_round(-5, 1), -3, "round(-2.5)=-3 (half away)");
    CHECK_EQ((long long)fxp_rshift_round(3, -2), 12, "negative shift = left shift");
    CHECK_EQ((long long)fxp_rshift_round(INT64_MIN, 1), -((long long)1 << 62), "rshift_round INT64_MIN");

    /* ---------------- conversion ---------------- */
    CHECK_EQ(fxp_from_float(1.0f), ONE, "from_float 1.0");
    CHECK_EQ(fxp_from_float(-1.5f), -3 * ONE / 2, "from_float -1.5");
    CHECK_EQ(fxp_from_float(1.0e9f), INT32_MAX, "from_float saturates");
    CHECK_EQ(fxp_from_float(-1.0e9f), INT32_MIN, "from_float saturates low");
    CHECK_NEAR(fxp_to_float(fxp_from_float(0.1f)), 0.1, 0.5 / ONE, "round trip 0.1 within 1/2 LSB");
    CHECK_EQ(fxp_from_double_q(0.0 / 0.0, 32, 16), 0, "NaN -> 0");
    CHECK_EQ(fxp_from_double_q(1.0, 16, 15), 32767, "+1.0 in Q1.15 saturates (no wrap-around)");
    CHECK_EQ(fxp_from_double_q(-1.0, 16, 15), -32768, "-1.0 in Q1.15 exact");
    CHECK_EQ(fxp_from_double_q(1.0, 16, 14), 16384, "+1.0 in Q2.14 exact");
    CHECK_EQ(fxp_from_double_q(1.0, 32, 30), 1 << 30, "+1.0 in Q2.30 exact");
    CHECK_EQ(fxp_from_double_q(0.5 / 65536, 32, 16), 1, "half LSB rounds away (+)");
    CHECK_EQ(fxp_from_double_q(-0.5 / 65536, 32, 16), -1, "half LSB rounds away (-)");

    /* ---------------- randomized property tests vs double ---------------- */
    fxp_stats_reset();
    for (i = 0; i < 100000; i++) {
        int32_t a = (int32_t)(t_rand() >> 12) - (1 << 19); /* about +-8 in Q16.16 */
        int32_t b = (int32_t)(t_rand() >> 12) - (1 << 19);
        double ex = ((double)a / ONE) * ((double)b / ONE);
        int32_t r = fxp_mul_round(a, b);
        int32_t f = fxp_mul(a, b);
        double e = fabs((double)r / ONE - ex);
        double ef = (double)f / ONE - ex;
        if (e > 0.5 / ONE + 1e-12) { CHECK(0, "mul_round error > 1/2 LSB"); break; }
        if (!(ef <= 1e-12 && ef > -1.0 / ONE - 1e-12)) { CHECK(0, "floor mul error outside (-1 LSB, 0]"); break; }
        if (b != 0) {
            double exd = ((double)a / ONE) / ((double)b / ONE);
            if (fabs(exd) < 32000.0) {
                int32_t d = fxp_div(a, b);
                if (fabs((double)d / ONE - exd) > 1.0 / ONE + 1e-12) { CHECK(0, "div error > 1 LSB"); break; }
            }
        }
    }
    CHECK(1, "100k randomized mul/div property checks done");
    return TEST_END();
}
