/*
 * s4_compare.h -- runs Model A (float) and Model B (fixed) in lock step on a
 * synthetic scenario and accumulates the numerical-accuracy metrics.
 * Analysis tooling (double allowed).  Not reentrant (static frame buffers).
 */
#ifndef S4_COMPARE_H
#define S4_COMPARE_H

#include <stdio.h>
#include "ekf_fixed_divider.h"
#include "ekf_fixed_math.h"
#include "ekf_s4_models.h"
#include "s4_dataset.h"
#include "s4_metrics.h"

#define S4_NGRP 14
#define S4_LOG_EVERY 50 /* log tick = 0.1 s at 500 Hz */

extern const char *const S4_GRP_NAME[S4_NGRP];
extern const char *const S4_GRP_UNIT[S4_NGRP];

typedef struct {
    s4_acc_t err[S4_NGRP];  /* fixed model vs float reference                */
    s4_acc_t qerr[S4_NGRP]; /* float state -> fixed -> back (pure quantization) */
    s4_acc_t tr_f[3];       /* float model vs TRUTH: position, velocity, attitude  */
    s4_acc_t tr_x[3];       /* fixed model vs TRUTH */
    s4_acc_t sens[3];       /* float model vs float model whose IMU inputs are perturbed by
                               ~1 float32 epsilon (relative): the reference's own sensitivity floor
                               for position, velocity, attitude */
    fxp_counters_t cnt[FXP_GRP_COUNT];
    uint32_t pd_clamp, bierman_abort, float_pd_viol;
    double max_drift_f, max_drift_x; /* max |p_i| over all frames [m] */
    double min_D_f, min_D_x;         /* min normalized D over log ticks  */
    double max_Un, max_Dn, peak_dx_n, max_abs_p, max_abs_v; /* range analysis (float reference, 0.1 s sampling) */
    double qnorm_dev_x;              /* max | |q_fixed| - 1 | over all frames (order-1 Taylor normalization) */
    int nan_flag;
    int n_frames, n_log;
} s4_cmp_t;

int s4_cmp_run(s4_scenario_t sc, double duration_s, int cfg_idx, const s4_params_t *prm,
               div_mode_t divmode, FILE *ts_csv, s4_cmp_t *out);

/* Same, with an explicit configuration (for design-space sweeps). */
int s4_cmp_run_cfg(s4_scenario_t sc, double duration_s, const ekf_qcfg_t *cfg,
                   const s4_params_t *prm, div_mode_t divmode, FILE *ts_csv, s4_cmp_t *out);

#endif /* S4_COMPARE_H */
