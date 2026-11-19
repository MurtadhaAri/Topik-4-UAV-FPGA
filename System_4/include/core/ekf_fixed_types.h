/*
 * ekf_fixed_types.h -- fixed-point mirrors of the ICD structures.
 *
 * The layout/indices mirror ekf_types.h exactly (NOM_DIM=16, ERR_DIM=15,
 * U_UPPER_SIZE=120, U_IDX unchanged).  Only the element type differs:
 * every element is an int32_t holding a signed value of at most `wl` bits.
 * The Q-format (fractional bits) of each field group is described by
 * ekf_qcfg_t and is NOT the same for all groups.
 *
 * Normalization: dx, U and D are stored in the NORMALIZED domain
 *   dx_n[i]  = dx[i] / s_i
 *   U_n[i,j] = U[i,j] * s_j / s_i        (unit diagonal preserved)
 *   D_n[i]   = D[i]   / s_i^2
 * with s_i = 2^s4_scale_log2[i] (see ekf_fixed_config.h).  This is a
 * similarity transform of P = U D U^T, so UDU^T structure is preserved.
 */
#ifndef EKF_FIXED_TYPES_H
#define EKF_FIXED_TYPES_H

#include <stdint.h>
#include "ekf_types.h"
#include "ekf_sensor_packet.h"
#include "ekf_fixed_config.h"

typedef struct {
    uint8_t wl;   /* word length in bits (<= 32), value is sign-extended in int32 */
    uint8_t frac; /* fractional bits */
} fxp_fmt_t;

typedef struct {
    int32_t p[3];
    int32_t v[3];
    int32_t q[4];
    int32_t ba[3];
    int32_t bg[3];
} nominal_state_fixed_t;

typedef struct {
    nominal_state_fixed_t x_nom;
    int32_t dx[ERR_DIM];     /* normalized error state                  */
    int32_t U[U_UPPER_SIZE]; /* normalized unit upper-triangular, U_IDX */
    int32_t D[ERR_DIM];      /* normalized diagonal                     */
} es_ekf_fixed_t;

/* One complete quantization configuration (PROPOSED design parameters). */
typedef struct {
    const char *name;
    /* nominal path (ICD Role 4: 32-bit fixed or float) */
    fxp_fmt_t p, v, q, b;
    /* error / covariance path (ICD Role 4: 16..24 bit after normalization) */
    fxp_fmt_t dx, U, D, Phi;
} ekf_qcfg_t;

/* Fixed-point sensor packet.  Formats:
 *   dt        Q2.30
 *   imu_acc   Q8.24 [m/s^2]      imu_gyr  Q8.24 [rad/s]
 *   meas[0..2] position/baro in cfg.p format, meas[3..5] velocity in cfg.v format
 */
typedef struct {
    uint32_t timestamp_us;
    int32_t dt;
    sensor_type_e type;
    int32_t imu_acc[3];
    int32_t imu_gyr[3];
    int32_t meas[6];
} sensor_frame_fixed_t;

/* Study configurations used by the analysis (see docs/Q_FORMAT_ANALYSIS.md) */
#define S4_NUM_CFG 3
const ekf_qcfg_t *s4_cfg_get(int index); /* 0:W16 1:W24 2:W32 */

#endif /* EKF_FIXED_TYPES_H */
