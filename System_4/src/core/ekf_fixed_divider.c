/*
 * ekf_fixed_divider.c -- baseline / restoring / non-restoring dividers.
 * Fixed loop bounds (nbits), no recursion, no dynamic memory, no FP.
 */
#include "ekf_fixed_divider.h"

#include "ekf_fixed_math.h"

static div_mode_t g_pipe_mode = DIV_BASELINE;

void div_set_pipeline_mode(div_mode_t mode) { g_pipe_mode = mode; }
div_mode_t div_get_pipeline_mode(void) { return g_pipe_mode; }

const char *div_mode_name(div_mode_t mode)
{
    switch (mode) {
    case DIV_BASELINE: return "baseline";
    case DIV_RESTORING: return "restoring";
    case DIV_NONRESTORING: return "non-restoring";
    default: return "?";
    }
}

uint8_t div_bit_width(uint64_t n)
{
    uint8_t w = 0;
    while (n != 0u) {
        w++;
        n >>= 1;
    }
    return w;
}

/* Restoring: shift one dividend bit into the partial remainder, trial
 * subtract, restore (keep old value) when the result would be negative.
 * Requires d != 0 and d <= 2^63 so that 2*r+1 cannot overflow uint64. */
uint64_t divu_restoring(uint64_t n, uint64_t d, uint8_t nbits, uint64_t *rem)
{
    uint64_t q = 0, r = 0;
    int i;
    if (nbits > 64u) nbits = 64u;
    for (i = (int)nbits - 1; i >= 0; i--) {
        r = (r << 1) | ((n >> i) & 1u);
        if (r >= d) {
            r -= d;
            q |= (uint64_t)1 << i;
        }
    }
    if (rem) *rem = r;
    return q;
}

/* Non-restoring: the partial remainder may go negative; the next step adds
 * the divisor instead of restoring.  One final correction if r < 0.
 * Partial remainder is kept as a signed value with |r| < 2*d (d < 2^62). */
uint64_t divu_nonrestoring(uint64_t n, uint64_t d, uint8_t nbits, uint64_t *rem)
{
    uint64_t q = 0;
    int64_t r = 0;
    int64_t sd = (int64_t)d;
    int i;
    if (nbits > 64u) nbits = 64u;
    for (i = (int)nbits - 1; i >= 0; i--) {
        int64_t bit = (int64_t)((n >> i) & 1u);
        if (r >= 0)
            r = 2 * r + bit - sd;
        else
            r = 2 * r + bit + sd;
        if (r >= 0) q |= (uint64_t)1 << i;
    }
    if (r < 0) r += sd; /* correction step */
    if (rem) *rem = (uint64_t)r;
    return q;
}

static uint64_t mag64(int64_t v) { return (v < 0) ? (uint64_t)0 - (uint64_t)v : (uint64_t)v; }

static uint64_t divu_dispatch(uint64_t n, uint64_t d, div_mode_t mode, div_info_t *info,
                              uint64_t *rem)
{
    uint64_t q;
    uint8_t nb = div_bit_width(n);
    if (nb == 0u) nb = 1u;
    if (mode == DIV_RESTORING) {
        q = divu_restoring(n, d, nb, rem);
        if (info) {
            info->iterations = nb;
            info->cycles_est = nb;
        }
    } else if (mode == DIV_NONRESTORING) {
        q = divu_nonrestoring(n, d, nb, rem);
        if (info) {
            info->iterations = nb;
            info->cycles_est = (uint32_t)nb + 1u;
        }
    } else {
        q = n / d;
        if (rem) *rem = n % d;
        if (info) {
            info->iterations = 0;
            info->cycles_est = 0; /* software: not a hardware cycle count */
        }
    }
    return q;
}

static int64_t div_special(int64_t n, int64_t d, div_info_t *info)
{
    (void)d;
    if (info) {
        info->div0 = 1;
        info->iterations = 0;
        info->cycles_est = 0;
    }
    fxp_stats_note_div0();
    if (n == 0) return 0;
    return (n > 0) ? INT64_MAX : -INT64_MAX;
}

int64_t divs64(int64_t n, int64_t d, div_mode_t mode, div_info_t *info)
{
    uint64_t q;
    int neg;
    if (info) info->div0 = 0;
    if (d == 0) return div_special(n, d, info);
    neg = ((n < 0) != (d < 0));
    q = divu_dispatch(mag64(n), mag64(d), mode, info, (uint64_t *)0);
    if (q > (uint64_t)INT64_MAX) q = (uint64_t)INT64_MAX;
    return neg ? -(int64_t)q : (int64_t)q;
}

int64_t divs64_round(int64_t n, int64_t d, div_mode_t mode, div_info_t *info)
{
    uint64_t q, r, md;
    int neg;
    if (info) info->div0 = 0;
    if (d == 0) return div_special(n, d, info);
    neg = ((n < 0) != (d < 0));
    md = mag64(d);
    q = divu_dispatch(mag64(n), md, mode, info, &r);
    if (r >= md - r) q++; /* 2r >= d  <=>  r >= d - r  (no overflow) */
    if (q > (uint64_t)INT64_MAX) q = (uint64_t)INT64_MAX;
    return neg ? -(int64_t)q : (int64_t)q;
}

int32_t div_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr,
              div_mode_t mode, div_info_t *info)
{
    int64_t num = (int64_t)a;
    int64_t den = (int64_t)b;
    int s = (int)fr + (int)fb - (int)fa;
    if (info) info->div0 = 0;
    if (b == 0) {
        (void)div_special(num, den, info);
        if (a == 0) return 0;
        fxp_stats_note_sat(0);
        return (a > 0) ? INT32_MAX : INT32_MIN;
    }
    if (s >= 0) {
        num = fxp_shl64_sat(num, s);
    } else if (-s > 32) {
        if (info) {
            info->iterations = 0;
            info->cycles_est = 0;
        }
        return 0;
    } else {
        den = fxp_shl64_sat(den, -s);
    }
    return fxp_saturate(divs64(num, den, mode, info));
}
