/*
 * ekf_fixed_divider.h -- division modules.
 *
 *   DIV_BASELINE       software divider (native int64 '/'), the reference.
 *   DIV_RESTORING      radix-2 restoring shift-and-subtract divider.
 *   DIV_NONRESTORING   radix-2 non-restoring divider (+1 correction step).
 *
 * NOT implemented: SRT (radix-4 / redundant quotient digits).  It remains a
 * candidate for later work and must not be claimed from these files.
 *
 * All three are bit-exact: floor(|n|/|d|) with the sign applied afterwards
 * (truncation toward zero).  They differ only in structure/iteration count.
 */
#ifndef EKF_FIXED_DIVIDER_H
#define EKF_FIXED_DIVIDER_H

#include <stdint.h>

typedef enum { DIV_BASELINE = 0, DIV_RESTORING = 1, DIV_NONRESTORING = 2 } div_mode_t;

typedef struct {
    uint32_t iterations; /* radix-2 iterations performed (= cycles, 1 bit/cycle) */
    uint32_t cycles_est; /* iterations (+1 correction step for non-restoring)    */
    int div0;            /* divisor was zero                                     */
} div_info_t;

/* Unsigned core: quotient = floor(n/d), nbits fixed iterations (<= 64). */
uint64_t divu_restoring(uint64_t n, uint64_t d, uint8_t nbits, uint64_t *rem);
uint64_t divu_nonrestoring(uint64_t n, uint64_t d, uint8_t nbits, uint64_t *rem);

/* Number of significant bits of n (0 for n==0): the fixed bound a hardware
 * divider of that dividend width would use. */
uint8_t div_bit_width(uint64_t n);

/* Signed 64/64 divide, truncation toward zero.  d==0 returns +/-INT64_MAX
 * (0 if n==0) and sets info->div0. */
int64_t divs64(int64_t n, int64_t d, div_mode_t mode, div_info_t *info);

/* Same but rounds half away from zero (used by the study pipeline). */
int64_t divs64_round(int64_t n, int64_t d, div_mode_t mode, div_info_t *info);

/* Selects the divider used by divs64_round()/fxp_div_wide() in the pipeline. */
void div_set_pipeline_mode(div_mode_t mode);
div_mode_t div_get_pipeline_mode(void);

/* Q-format divide with the selected structure; value = (a/2^fa)/(b/2^fb) in
 * Qfr, saturated to int32 (counts in fxp stats).  For DIV_BASELINE this is
 * identical to fxp_div_q(). */
int32_t div_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr,
              div_mode_t mode, div_info_t *info);

const char *div_mode_name(div_mode_t mode);

#endif /* EKF_FIXED_DIVIDER_H */
