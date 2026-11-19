/*
 * ekf_s4_models.h -- System 4 STUDY PIPELINES (Model A float, Model B fixed).
 *
 * IMPORTANT SCOPE NOTE
 *   These are NOT the team's ES-EKF.  They are minimal pipelines that follow
 *   the equations written in the ICD (nominal propagation, Phi, Thornton MWGS,
 *   Bierman scalar update with unit-row H, error injection/reset) so that
 *   quantization can be measured on U, D, dx and the nominal state BEFORE the
 *   members 1-3 modules are available.  Status: PENDING INTEGRATION; when the
 *   real modules exist, replace these call-by-call (see INTERFACE_CONTRACT.md).
 *
 *   Supported measurements: GNSS pos+vel (6 scalar channels) and barometer.
 *   Magnetometer (SENS_MAG) is NOT implemented (PENDING).
 *
 *   Every process/measurement noise value below is a PROPOSED parameter
 *   except R and D0, which are the ICD baseline numbers.
 */
#ifndef EKF_S4_MODELS_H
#define EKF_S4_MODELS_H

#include "ekf_fixed_types.h"

/* State groups: 0 position, 1 velocity, 2 attitude, 3 accel bias, 4 gyro bias */
typedef struct {
    float D0[5];  /* ICD: 1.0, 0.5, 0.05, 0.01, 0.0001                         */
    float qc[5];  /* PROPOSED continuous noise density per group [unit^2/s]    */
    float r_pos;  /* ICD 0.25 m^2   */
    float r_vel;  /* ICD 0.04 (m/s)^2 */
    float r_baro; /* ICD 0.16 m^2   */
    int cov_decim;/* covariance time update every N IMU samples (PROPOSED)     */
} s4_params_t;

void s4_params_default(s4_params_t *p);

/* ============================ Model A: float32 reference ================= */
typedef struct {
    es_ekf_t ekf;
    s4_params_t prm;
    float dt_acc;
    int cnt;
    float a_unb[3], w_unb[3], R[9];
    uint32_t n_pd_viol; /* Thornton produced D<=0 (should stay 0) */
    float peak_dx_n;    /* max |dx_i|/s_i seen just before error injection (range analysis) */
} s4f_model_t;

void s4f_init(s4f_model_t *m, const s4_params_t *prm);
void s4f_step(s4f_model_t *m, const sensor_frame_t *f);
void s4f_set_velocity(s4f_model_t *m, const float v[3]); /* initial nominal velocity */

/* ============================ Model B: fixed-point ======================= */
typedef struct {
    es_ekf_fixed_t ekf;
    ekf_qcfg_t cfg;
    int decim;
    int32_t dt_acc; /* Q2.30 */
    int cnt;
    int32_t a_unb[3], w_unb[3]; /* Q8.24 */
    int32_t Rm[9];              /* Q2.30 */
    int64_t qc_n[5];            /* qc/s^2 per second, frac 50 */
    int32_t r_n[3];             /* R/s^2 in Q3.29: pos, vel, baro */
    int32_t g_q;                /* gravity Q8.24 */
    uint32_t n_pd_clamp;        /* D<=0 clamped to 1 LSB (loss of positive definiteness) */
    uint32_t n_bierman_abort;   /* alpha <= 0 */
} s4x_model_t;

void s4x_init(s4x_model_t *m, const ekf_qcfg_t *cfg, const s4_params_t *prm);
void s4x_step(s4x_model_t *m, const sensor_frame_fixed_t *f);
void s4x_set_velocity(s4x_model_t *m, const float v[3]); /* initial nominal velocity (I/O conversion) */

#endif /* EKF_S4_MODELS_H */
