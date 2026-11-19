/* CORDIC sin/cos against libm (reference only inside the test). */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "ekf_fixed_math.h"
#include "ekf_fixed_trig.h"
#include "test_util.h"

static double ang(double a) { return fxp_to_double_q(cordic_angle_from_double(a), 29); }

int main(void)
{
    const int N = 40001;
    int k;
    double ss = 0, sc = 0, ms = 0, mc = 0, mnorm = 0;
    int32_t s, c;
    TEST_BEGIN("ekf_fixed_trig (CORDIC)");

    CHECK_EQ(CORDIC_ITERATIONS, 30, "default iteration count");
    CHECK_EQ(CORDIC_ATAN_Q29[0], 421657428, "atan(1) table entry = pi/4 in Q3.29");
    CHECK_NEAR(fxp_to_double_q(CORDIC_K_Q30, 30), 0.6072529350088813, 1e-9, "gain constant K");
    CHECK_NEAR(fxp_to_double_q(CORDIC_PI_Q29, 29), 3.141592653589793, 1e-9, "pi constant");

    /* sweep over the full input range */
    for (k = 0; k < N; k++) {
        double a = -3.141592653589793 + 2.0 * 3.141592653589793 * k / (N - 1);
        double aq = ang(a), es, ec, nrm;
        cordic_sincos(cordic_angle_from_double(a), &s, &c);
        es = fabs(fxp_to_double_q(s, 30) - sin(aq));
        ec = fabs(fxp_to_double_q(c, 30) - cos(aq));
        nrm = fabs(fxp_to_double_q(s, 30) * fxp_to_double_q(s, 30) +
                   fxp_to_double_q(c, 30) * fxp_to_double_q(c, 30) - 1.0);
        ss += es * es;
        sc += ec * ec;
        if (es > ms) ms = es;
        if (ec > mc) mc = ec;
        if (nrm > mnorm) mnorm = nrm;
    }
    printf("  30 iterations: RMS sin %.3e cos %.3e | max sin %.3e cos %.3e | max |s^2+c^2-1| %.3e\n",
           sqrt(ss / N), sqrt(sc / N), ms, mc, mnorm);
    CHECK(sqrt(ss / N) < 1e-8 && sqrt(sc / N) < 1e-8, "RMS error < 1e-8 (about 10 LSB of Q2.30)");
    CHECK(ms < 4e-8 && mc < 4e-8, "max error < 4e-8 (about 43 LSB of Q2.30)");
    CHECK(mnorm < 1e-7, "unit-norm error < 1e-7");

    /* special angles, no wrap-around of +1 */
    cordic_sincos(0, &s, &c);
    CHECK(c > 0 && llabs((long long)c - (1 << 30)) < 64, "cos(0) = +1 without wrap-around");
    CHECK(llabs(s) < 64, "sin(0) ~ 0");
    cordic_sincos((int32_t)CORDIC_HALFPI_Q29, &s, &c);
    CHECK(s > 0 && llabs((long long)s - (1 << 30)) < 64, "sin(pi/2) = +1");
    CHECK(llabs(c) < 64, "cos(pi/2) ~ 0");
    cordic_sincos(-(int32_t)CORDIC_HALFPI_Q29, &s, &c);
    CHECK(s < 0 && llabs((long long)s + (1 << 30)) < 64, "sin(-pi/2) = -1");
    cordic_sincos((int32_t)CORDIC_PI_Q29, &s, &c);
    CHECK(c < 0 && llabs((long long)c + (1 << 30)) < 64, "cos(pi) = -1");
    CHECK(llabs(s) < 64, "sin(pi) ~ 0");
    cordic_sincos(-(int32_t)CORDIC_PI_Q29, &s, &c);
    CHECK(c < 0 && llabs((long long)c + (1 << 30)) < 64, "cos(-pi) = -1");

    /* quadrant signs */
    cordic_sincos(cordic_angle_from_double(2.0), &s, &c);
    CHECK(s > 0 && c < 0, "quadrant II");
    cordic_sincos(cordic_angle_from_double(-2.0), &s, &c);
    CHECK(s < 0 && c < 0, "quadrant III");
    cordic_sincos(cordic_angle_from_double(-1.0), &s, &c);
    CHECK(s < 0 && c > 0, "quadrant IV");

    /* inputs in (pi, 4) are wrapped by -2*pi exactly once (Q3.29 spans [-4,4)) */
    for (k = 0; k < 200; k++) {
        double a = 3.15 + (3.99 - 3.15) * k / 199.0, aq = ang(a);
        cordic_sincos(cordic_angle_from_double(a), &s, &c);
        if (fabs(fxp_to_double_q(s, 30) - sin(aq)) > 4e-8 || fabs(fxp_to_double_q(c, 30) - cos(aq)) > 4e-8) {
            CHECK(0, "wrapped angle in (pi,4) wrong");
            break;
        }
        a = -a;
        aq = ang(a);
        cordic_sincos(cordic_angle_from_double(a), &s, &c);
        if (fabs(fxp_to_double_q(s, 30) - sin(aq)) > 4e-8 || fabs(fxp_to_double_q(c, 30) - cos(aq)) > 4e-8) {
            CHECK(0, "wrapped angle in [-4,-pi) wrong");
            break;
        }
    }
    CHECK(1, "wrapped angle ranges checked");
    cordic_sincos(2147483647, &s, &c); /* ~ +3.99999 rad */
    CHECK_NEAR(fxp_to_double_q(s, 30), sin(3.9999999981), 4e-8, "max angle input sin");

    /* error shrinks with iteration count (about 2 bits per iteration) */
    {
        double e_prev = 1e9;
        int its[5] = {8, 12, 16, 20, 24}, i;
        for (i = 0; i < 5; i++) {
            double e = 0;
            for (k = 0; k < 2001; k++) {
                double a = -3.0 + 6.0 * k / 2000.0, aq = ang(a);
                cordic_sincos_n(cordic_angle_from_double(a), (uint8_t)its[i], &s, &c);
                if (fabs(fxp_to_double_q(s, 30) - sin(aq)) > e) e = fabs(fxp_to_double_q(s, 30) - sin(aq));
            }
            printf("  N=%2d  max sin error %.3e\n", its[i], e);
            CHECK(e < e_prev, "error decreases with more iterations");
            e_prev = e;
        }
    }
    /* iteration request above the maximum is clamped, not an out-of-bounds table read */
    cordic_sincos_n(cordic_angle_from_double(1.0), 200, &s, &c);
    CHECK_NEAR(fxp_to_double_q(s, 30), sin(1.0), 4e-8, "N clamped to maximum");
    return TEST_END();
}
