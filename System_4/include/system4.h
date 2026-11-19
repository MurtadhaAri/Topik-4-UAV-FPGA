/*
 * system4.h -- single entry header of System 4
 *              "Analis Kuantisasi Fixed-Point & Bit-Accurate Model"
 *
 * One include gives a team member:
 *   (1) the BEST quantization configuration (s4_best_cfg()),
 *   (2) the acceptance thresholds that configuration must meet,
 *   (3) the DUT (device-under-test) template: the 3-function interface a
 *       fixed-point ES-EKF implementation must provide so System 4's proof
 *       harness can compare it against the float32 reference, plus the
 *       canonical stage order of one filter step.
 *
 * Status labels used below:  ICD = agreed in the team documents,
 *                            PROPOSED = System 4 design parameter (not agreed),
 *                            MEASURED = produced by tests/proof/system4_proof.c.
 * Notation: Qm.n = m integer bits INCLUDING sign, n fractional bits.
 */
#ifndef SYSTEM4_H
#define SYSTEM4_H

#include <stddef.h>
#include <stdint.h>

#include "ekf_fixed_adapter.h"
#include "ekf_fixed_divider.h"
#include "ekf_fixed_math.h"
#include "ekf_fixed_trig.h"
#include "ekf_fixed_types.h"
#include "ekf_s4_models.h"

/* ======================================================================== *
 * 1. BEST CONFIGURATION  (PROPOSED, selected and checked by system4_proof)
 *
 *   nominal path (32-bit)      p Q16.16   v Q8.24   q Q2.30   bias Q4.28   (ICD Role 4 + System 4)
 *   error path   (storage)     dx Q4.20 (24 b)   U Q7.17 (24 b)
 *                              D  Q1.31 (32 b)   dPhi Q1.23 (24 b)
 *   working registers          W Q8.24 (int32), alpha/g/k Q3.29, int64 accumulators
 *   normalization scale s_i    16,16,16, 16,16,16, 0.5,0.5,0.5, 2,2,2, 0.25,0.25,0.25
 *   quaternion normalization   Newton inverse-sqrt (NOT the order-1 Taylor rule)
 *
 *   Why: D storage width is the dominant accuracy lever (ablation, docs/Q_FORMAT_ANALYSIS.md);
 *   dx, dPhi at 24 bit are indistinguishable from 32 bit; U at 24 bit costs attitude accuracy.
 * ======================================================================== */
static inline const ekf_qcfg_t *s4_best_cfg(void)
{
    static const ekf_qcfg_t best = {
        "BEST",
        {32, 16}, {32, 24}, {32, 30}, {32, 28}, /* p, v, q, bias            */
        {24, 20}, {24, 17}, {32, 31}, {24, 23}  /* dx, U, D, dPhi           */
    };
    return &best;
}

/* Storage cost of the error path in bits: dx[15] + U[120] + D[15]. */
static inline uint32_t s4_cfg_storage_bits(const ekf_qcfg_t *c)
{
    return 15u * c->dx.wl + 120u * c->U.wl + 15u * c->D.wl;
}

/* ======================================================================== *
 * 2. ACCEPTANCE THRESHOLDS (PROPOSED -- the ICD only defines the static
 *    criterion).  Fixed model vs float32 reference, RMSE over a 60 s run.
 * ======================================================================== */
#define S4_ACC_STATIC_DRIFT_MM 1.0   /* ICD Role 5: position drift <= 1 mm      */
#define S4_ACC_POS_RMSE_M      5e-3  /* constvel, dynamic, highrange            */
#define S4_ACC_VEL_RMSE_MPS    5e-3
#define S4_ACC_ATT_RMSE_RAD    5e-3
#define S4_ACC_QNORM_DEV       1e-5  /* max | |q| - 1 |                         */
#define S4_ACC_MAX_SATURATIONS 0     /* no saturation / overflow allowed        */
#define S4_ACC_MAX_PD_CLAMPS   0     /* D must stay positive (ICD: D_i > 0)     */

/* ======================================================================== *
 * 3. DUT TEMPLATE -- plug a fixed-point ES-EKF into the System 4 proof
 *
 *  A DUT is three functions + the size of its private context (static memory,
 *  no malloc).  Formats follow ekf_fixed_types.h: es_ekf_fixed_t holds the
 *  NORMALIZED dx/U/D and the nominal state in the cfg formats; sensor_frame_fixed_t
 *  carries dt Q2.30, imu Q8.24, meas in p/v formats.
 *
 *  Canonical step (every DUT must keep this order -- ICD Role 1-3):
 *     1. nominal propagation        a_unb = a_m - b_a ; p,v,q update ; q normalised
 *     2. covariance time update     dPhi from R, a_unb, w_unb ; Thornton MWGS (U, D)
 *     3. per scalar channel         innovation y = z - h(x_nom) ; bierman(i0, y, R)
 *     4. error injection + reset    x_nom += dx ; q = norm(q (x) [1, dtheta/2]) ; dx = 0
 *  Rules: static arrays, flattened 1-D, fixed loop bounds, no sqrt/inverse in the
 *  covariance loop, saturating arithmetic from ekf_fixed_math.h only.
 *  The study implementation (S4_DUT_STUDY) shows one way; S4_DUT_TEMPLATE is a
 *  compilable skeleton with TODO stages that a member fills in.
 * ======================================================================== */
#define S4_DUT_CTX_MAX 65536u

typedef struct {
    const char *name;
    size_t ctx_size; /* must be <= S4_DUT_CTX_MAX */
    void (*init)(void *ctx, const ekf_qcfg_t *cfg, const s4_params_t *prm, const float v0[3]);
    void (*step)(void *ctx, const sensor_frame_fixed_t *frame);
    void (*export_state)(const void *ctx, es_ekf_fixed_t *out);
    /* optional diagnostics: number of "D <= 0 clamped" events (may be NULL) */
    uint32_t (*pd_clamps)(const void *ctx);
} s4_dut_t;

extern const s4_dut_t S4_DUT_STUDY;    /* System 4 study pipeline (reference DUT)        */
extern const s4_dut_t S4_DUT_TEMPLATE; /* skeleton: stages are TODO stubs -> must FAIL   */

#endif /* SYSTEM4_H */
