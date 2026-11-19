/*
 * ekf_fixed_adapter.h -- conversion between the ICD float structures and the
 * fixed-point mirrors (including (de)normalization).
 *
 * Deviation from the task brief: the converters take an extra `cfg` argument
 * because the Q-format is configurable (Q16.16, Q2.30, 16/24/32 bit ...).
 * Conversion uses float/double and is therefore ONLY for I/O and analysis,
 * never for the fixed-point core.
 */
#ifndef EKF_FIXED_ADAPTER_H
#define EKF_FIXED_ADAPTER_H

#include "ekf_fixed_types.h"

void ekf_float_to_fixed(const es_ekf_t *src, es_ekf_fixed_t *dst, const ekf_qcfg_t *cfg);
void ekf_fixed_to_float(const es_ekf_fixed_t *src, es_ekf_t *dst, const ekf_qcfg_t *cfg);

void sensor_frame_to_fixed(const sensor_frame_t *src, sensor_frame_fixed_t *dst,
                           const ekf_qcfg_t *cfg);

/* Scale factor s_i = 2^s4_scale_log2[i] as double (analysis helper). */
double s4_scale(int i);

/* Normalized covariance entries as double (for analysis / comparison):
 *   U_n(i,j) = U(i,j)*s_j/s_i ,  D_n(i) = D(i)/s_i^2                       */
double s4_float_Un(const es_ekf_t *e, int i, int j);
double s4_float_Dn(const es_ekf_t *e, int i);
double s4_fixed_Un(const es_ekf_fixed_t *e, const ekf_qcfg_t *cfg, int i, int j);
double s4_fixed_Dn(const es_ekf_fixed_t *e, const ekf_qcfg_t *cfg, int i);

#endif /* EKF_FIXED_ADAPTER_H */
