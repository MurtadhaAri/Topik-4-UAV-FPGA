/*
 * ekf_fixed_math.h -- bit-accurate saturating fixed-point arithmetic.
 *
 * Containers are int32_t; every product uses an int64_t intermediate.
 * "wl" arguments model hardware words narrower than the container (16/24/32).
 * No floating point is used by the arithmetic functions; the only float/double
 * users are the explicit conversion helpers at the bottom.
 *
 * Rounding conventions (documented, tested):
 *   fxp_mul*, fxp_shift_right ... "floor"  (arithmetic shift, rounds to -inf)
 *   *_round                   ... round half away from zero (symmetric)
 *   fxp_div                   ... truncation toward zero (C99 '/')
 */
#ifndef EKF_FIXED_MATH_H
#define EKF_FIXED_MATH_H

#include <stdint.h>
#include "ekf_fixed_config.h"

/* ------------------------------------------------------------------ stats */
typedef enum {
    FXP_GRP_MISC = 0,
    FXP_GRP_NOM_P,
    FXP_GRP_NOM_V,
    FXP_GRP_NOM_Q,
    FXP_GRP_NOM_B,
    FXP_GRP_DX,
    FXP_GRP_U,
    FXP_GRP_D,
    FXP_GRP_PHI,
    FXP_GRP_WORK,
    FXP_GRP_COUNT
} fxp_group_t;

typedef struct {
    uint32_t sat;       /* result clipped to the target word length          */
    uint32_t ovf;       /* pre-clip value did not even fit int32 (container) */
    uint32_t underflow; /* non-zero operands produced an exact 0 result      */
    uint32_t div0;      /* division by zero                                  */
} fxp_counters_t;

void fxp_stats_reset(void);
void fxp_stats_set_group(fxp_group_t g);
fxp_counters_t fxp_stats_get(fxp_group_t g);
fxp_counters_t fxp_stats_total(void);
const char *fxp_group_name(fxp_group_t g);
void fxp_stats_note_div0(void);
void fxp_stats_note_sat(int container_overflow);
void fxp_stats_note_underflow(void);

/* ------------------------------------------------------- format utilities */
int64_t fxp_wl_max(uint8_t wl); /*  2^(wl-1)-1 */
int64_t fxp_wl_min(uint8_t wl); /* -2^(wl-1)   */

/* ------------------------------------------------- 64-bit helper primitives */
int64_t fxp_asr64(int64_t v, int sh);       /* floor(v / 2^sh), sh >= 0        */
int64_t fxp_shl64_sat(int64_t v, int sh);   /* v * 2^sh saturated to int64     */
int64_t fxp_rshift_round(int64_t v, int sh);/* round-half-away; sh<0 = left shift */
int64_t fxp_add64_sat(int64_t a, int64_t b);
int64_t fxp_mul64_sat(int64_t a, int64_t b);

/* ------------------------------------------------ required public interface */
int32_t fxp_saturate(int64_t value);               /* -> int32, counted       */
int32_t fxp_sat_wl(int64_t value, uint8_t wl);     /* -> signed wl-bit, counted */
int32_t fxp_add_sat(int32_t a, int32_t b);
int32_t fxp_sub_sat(int32_t a, int32_t b);
int32_t fxp_mul(int32_t a, int32_t b);             /* default frac, floor     */
int32_t fxp_mul_round(int32_t a, int32_t b);       /* default frac, rounded   */
int32_t fxp_div(int32_t a, int32_t b);             /* default frac, trunc     */
int32_t fxp_shift_right(int32_t value, uint8_t shift); /* arithmetic, floor   */
int32_t fxp_from_float(float value);               /* default frac, rounded   */
float fxp_to_float(int32_t value);

#define MUL_FXP(a, b) fxp_mul((a), (b))

/* ------------------------- format-aware versions (fa, fb, fr = frac bits) */
int32_t fxp_mul_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr);
int32_t fxp_mul_round_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr);
int32_t fxp_div_q(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr);
int32_t fxp_mul_round_wl(int32_t a, int32_t b, uint8_t fa, uint8_t fb, uint8_t fr, uint8_t wl);
int32_t fxp_add_sat_wl(int32_t a, int32_t b, uint8_t wl);

/* ------------------------------------------------ conversion (I/O only) */
int32_t fxp_from_double_q(double v, uint8_t wl, uint8_t frac);
double fxp_to_double_q(int32_t raw, uint8_t frac);

#endif /* EKF_FIXED_MATH_H */
