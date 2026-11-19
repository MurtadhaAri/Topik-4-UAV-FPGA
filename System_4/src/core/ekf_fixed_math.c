/*
 * ekf_fixed_math.c -- bit-accurate saturating fixed-point arithmetic.
 * See ekf_fixed_math.h for conventions.  No malloc, no recursion, no FP in
 * the arithmetic routines.
 */
#include "ekf_fixed_math.h"

#include <math.h> /* conversion helpers only */
#include <string.h>

/* ------------------------------------------------------------ statistics */
static fxp_counters_t g_cnt[FXP_GRP_COUNT];
static fxp_group_t g_cur = FXP_GRP_MISC;

static const char *const k_group_names[FXP_GRP_COUNT] = {
    "misc", "nom_p", "nom_v", "nom_q", "nom_b", "dx", "U", "D", "Phi", "work"};

void fxp_stats_reset(void)
{
    memset(g_cnt, 0, sizeof g_cnt);
    g_cur = FXP_GRP_MISC;
}

void fxp_stats_set_group(fxp_group_t g)
{
    g_cur = ((int)g >= 0 && g < FXP_GRP_COUNT) ? g : FXP_GRP_MISC;
}

fxp_counters_t fxp_stats_get(fxp_group_t g)
{
    fxp_counters_t z = {0, 0, 0, 0};
    return ((int)g >= 0 && g < FXP_GRP_COUNT) ? g_cnt[g] : z;
}

fxp_counters_t fxp_stats_total(void)
{
    fxp_counters_t t = {0, 0, 0, 0};
    int i;
    for (i = 0; i < (int)FXP_GRP_COUNT; i++) {
        t.sat += g_cnt[i].sat;
        t.ovf += g_cnt[i].ovf;
        t.underflow += g_cnt[i].underflow;
        t.div0 += g_cnt[i].div0;
    }
    return t;
}

const char *fxp_group_name(fxp_group_t g)
{
    return ((int)g >= 0 && g < FXP_GRP_COUNT) ? k_group_names[g] : "?";
}

void fxp_stats_note_div0(void) { g_cnt[g_cur].div0++; }

void fxp_stats_note_underflow(void) { g_cnt[g_cur].underflow++; }

void fxp_stats_note_sat(int container_overflow)
{
    g_cnt[g_cur].sat++;
    if (container_overflow) g_cnt[g_cur].ovf++;
}

/* --------------------------------------------------------------- formats */
int64_t fxp_wl_max(uint8_t wl)
{
    if (wl < 2u) wl = 2u;
    if (wl > 63u) wl = 63u;
    return ((int64_t)1 << (wl - 1u)) - 1;
}

int64_t fxp_wl_min(uint8_t wl)
{
    if (wl < 2u) wl = 2u;
    if (wl > 63u) wl = 63u;
    return -((int64_t)1 << (wl - 1u));
}

/* ------------------------------------------------- 64-bit primitives */
/* Portable arithmetic shift (floor division by 2^sh), no reliance on
 * implementation-defined >> of negative values. */
int64_t fxp_asr64(int64_t v, int sh)
{
    if (sh <= 0) return v;
    if (sh >= 63) return (v < 0) ? -1 : 0;
    if (v >= 0) return v >> sh;
    return ~((~v) >> sh);
}

int64_t fxp_shl64_sat(int64_t v, int sh)
{
    int64_t lim;
    if (sh <= 0 || v == 0) return v;
    if (sh >= 63) return (v > 0) ? INT64_MAX : INT64_MIN;
    lim = INT64_MAX >> sh;
    if (v > lim) return INT64_MAX;
    if (v < -lim) return INT64_MIN;
    return (int64_t)((uint64_t)v << sh);
}

/* Round half away from zero. sh < 0 means a (saturating) left shift. */
int64_t fxp_rshift_round(int64_t v, int sh)
{
    uint64_t mag, half, res;
    if (sh <= 0) return fxp_shl64_sat(v, -sh);
    if (sh >= 63) return 0;
    mag = (v < 0) ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    half = (uint64_t)1 << (sh - 1);
    res = (mag + half) >> sh;
    return (v < 0) ? -(int64_t)res : (int64_t)res;
}

int64_t fxp_add64_sat(int64_t a, int64_t b)
{
    if (b > 0 && a > INT64_MAX - b) return INT64_MAX;
    if (b < 0 && a < INT64_MIN - b) return INT64_MIN;
    return a + b;
}

int64_t fxp_mul64_sat(int64_t a, int64_t b)
{
    uint64_t ma, mb;
    int neg;
    if (a == 0 || b == 0) return 0;
    neg = ((a < 0) != (b < 0));
    ma = (a < 0) ? (uint64_t)0 - (uint64_t)a : (uint64_t)a;
    mb = (b < 0) ? (uint64_t)0 - (uint64_t)b : (uint64_t)b;
    if (ma > (uint64_t)INT64_MAX / mb) return neg ? INT64_MIN : INT64_MAX;
    return neg ? -(int64_t)(ma * mb) : (int64_t)(ma * mb);
}

/* ------------------------------------------------------------ saturation */
int32_t fxp_sat_wl(int64_t v, uint8_t wl)
{
    int64_t hi = fxp_wl_max(wl);
    int64_t lo = fxp_wl_min(wl);
    if (v > hi) {
        fxp_stats_note_sat(v > (int64_t)INT32_MAX);
        return (int32_t)hi;
    }
    if (v < lo) {
        fxp_stats_note_sat(v < (int64_t)INT32_MIN);
        return (int32_t)lo;
    }
    return (int32_t)v;
}

int32_t fxp_saturate(int64_t v) { return fxp_sat_wl(v, 32u); }

int32_t fxp_add_sat(int32_t a, int32_t b) { return fxp_saturate((int64_t)a + (int64_t)b); }
int32_t fxp_sub_sat(int32_t a, int32_t b) { return fxp_saturate((int64_t)a - (int64_t)b); }
int32_t fxp_add_sat_wl(int32_t a, int32_t b, uint8_t wl)
{
    return fxp_sat_wl((int64_t)a + (int64_t)b, wl);
}

int32_t fxp_shift_right(int32_t value, uint8_t shift)
{
    return (int32_t)fxp_asr64((int64_t)value, (int)shift);
}

/* -------------------------------------------------------- multiplication */
int32_t fxp_mul_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr)
{
    int64_t prod = (int64_t)a * (int64_t)b; /* |prod| <= 2^62: never overflows */
    int sh = (int)fa + (int)fb - (int)fr;
    int64_t r = (sh >= 0) ? fxp_asr64(prod, sh) : fxp_shl64_sat(prod, -sh);
    if (r == 0 && a != 0 && b != 0) g_cnt[g_cur].underflow++;
    return fxp_saturate(r);
}

int32_t fxp_mul_round_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr)
{
    return fxp_mul_round_wl(a, b, fa, fb, fr, 32u);
}

int32_t fxp_mul_round_wl(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr, uint8_t wl)
{
    int64_t prod = (int64_t)a * (int64_t)b;
    int sh = (int)fa + (int)fb - (int)fr;
    int64_t r = fxp_rshift_round(prod, sh);
    if (r == 0 && a != 0 && b != 0) g_cnt[g_cur].underflow++;
    return fxp_sat_wl(r, wl);
}

int32_t fxp_mul(int32_t a, int32_t b)
{
    return fxp_mul_q(a, b, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC);
}

int32_t fxp_mul_round(int32_t a, int32_t b)
{
    return fxp_mul_round_q(a, b, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC);
}

/* -------------------------------------------------- baseline division */
/* value = (a/2^fa) / (b/2^fb) expressed in Qfr:  raw = a * 2^(fr+fb-fa) / b  */
int32_t fxp_div_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr)
{
    int64_t num = (int64_t)a;
    int64_t den = (int64_t)b;
    int s = (int)fr + (int)fb - (int)fa;
    if (b == 0) {
        fxp_stats_note_div0();
        if (a == 0) return 0;
        fxp_stats_note_sat(0);
        return (a > 0) ? INT32_MAX : INT32_MIN;
    }
    if (s >= 0) {
        num = fxp_shl64_sat(num, s);
    } else if (-s > 32) {
        return 0; /* |num| < 2^31 <= |den| * 2^(-s): quotient truncates to 0 */
    } else {
        den = fxp_shl64_sat(den, -s); /* |den| < 2^31, shift <= 32: fits */
    }
    return fxp_saturate(num / den); /* C99: truncation toward zero */
}

int32_t fxp_div(int32_t a, int32_t b)
{
    return fxp_div_q(a, b, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC, FXP_DEFAULT_FRAC);
}

/* ------------------------------------------------ conversions (I/O only) */
/* BEGIN-CONVERSION */
int32_t fxp_from_double_q(double v, uint8_t wl, uint8_t frac)
{
    double x;
    if (v != v) { /* NaN */
        fxp_stats_note_sat(1);
        return 0;
    }
    x = ldexp(v, (int)frac);
    if (x >= 9.0e18) return fxp_sat_wl(INT64_MAX, wl);
    if (x <= -9.0e18) return fxp_sat_wl(INT64_MIN, wl);
    x = (x >= 0.0) ? floor(x + 0.5) : -floor(-x + 0.5); /* half away from zero */
    return fxp_sat_wl((int64_t)x, wl);
}

double fxp_to_double_q(int32_t raw, uint8_t frac) { return ldexp((double)raw, -(int)frac); }

int32_t fxp_from_float(float value)
{
    return fxp_from_double_q((double)value, 32u, FXP_DEFAULT_FRAC);
}

float fxp_to_float(int32_t value) { return (float)fxp_to_double_q(value, FXP_DEFAULT_FRAC); }
/* END-CONVERSION */
