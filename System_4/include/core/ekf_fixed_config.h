/*
 * ekf_fixed_config.h -- compile-time configuration of System 4.
 *
 * Notation used everywhere in System 4:  Qm.n  = signed two's complement,
 * m integer bits INCLUDING the sign bit, n fractional bits, word length m+n.
 * (Q1.15 = 16 bit, Q1.23 = 24 bit, Q2.30 = 32 bit, Q16.16 = 32 bit.)
 *
 * Everything marked PROPOSED is a System 4 design parameter and NOT part of
 * the team agreement (ICD).  Items marked ICD are taken from the documents.
 */
#ifndef EKF_FIXED_CONFIG_H
#define EKF_FIXED_CONFIG_H

#include <stdint.h>
#include "ekf_types.h" /* ICD: NOM_DIM, ERR_DIM, U_UPPER_SIZE, U_IDX */

/* ICD 1.3: g = 9.80665 m/s^2 (the System 4 task brief said 9.81; ICD wins). */
#define S4_GRAVITY_MPS2 9.80665

/* Default fractional bits used by the 2-operand API fxp_mul/fxp_div/MUL_FXP.
 * PROPOSED: Q16.16.  Override with -DFXP_DEFAULT_FRAC=n.                    */
#ifndef FXP_DEFAULT_FRAC
#define FXP_DEFAULT_FRAC 16u
#endif

/* CORDIC (PROPOSED): angle Q3.29 rad, outputs Q2.30, fixed iteration count. */
#ifndef CORDIC_ITERATIONS
#define CORDIC_ITERATIONS 30u
#endif
#define CORDIC_MAX_ITERATIONS 30u
#define CORDIC_ANGLE_FRAC 29u
#define CORDIC_TRIG_FRAC 30u

/* Internal fixed formats of the study pipeline (PROPOSED) */
#define S4_IMU_FRAC 24u   /* acc/gyr input and a_unb, w_unb: Q8.24 */
#define S4_DT_FRAC 30u    /* dt, rotation matrix R, quaternion: Q2.30 */
#define S4_WORK_FRAC 24u  /* W matrix of Thornton, innovations: Q8.24 */
#define S4_ALPHA_FRAC 29u /* Bierman alpha, g, k: Q3.29 */

/* Covariance time update every N IMU samples (PROPOSED, N=10 -> 50 Hz @500Hz) */
#define S4_COV_DECIM_DEFAULT 10

/* Normalization scale factors s_i = 2^S4_SCALE_LOG2[i]  (PROPOSED).
 * Chosen as the smallest power of two >= the ICD error bounds
 * (10 m, 10 m/s, 0.5 rad, 2 m/s^2, 0.2 rad/s) -> 16, 16, 0.5, 2, 0.25.
 * Power-of-two scales make (de)normalization a pure bit shift.             */
extern const int8_t s4_scale_log2[ERR_DIM];

#endif /* EKF_FIXED_CONFIG_H */
