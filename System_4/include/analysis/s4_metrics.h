/*
 * s4_metrics.h -- error accumulators used by the RMSE analysis.
 *   RMSE = sqrt(mean((ref - test)^2)), MAE, max |err|, mean relative error.
 * Analysis tooling (uses double); never linked into the fixed-point core.
 */
#ifndef S4_METRICS_H
#define S4_METRICS_H

#include <stdint.h>

typedef struct {
    uint64_t n;
    double sum_sq, sum_abs, max_abs;
    double sum_rel;
    uint64_t n_rel;
} s4_acc_t;

#define S4_REL_EPS 1e-9 /* |ref| below this is excluded from the relative error */

void s4_acc_init(s4_acc_t *a);
void s4_acc_add(s4_acc_t *a, double ref, double test);
double s4_acc_rmse(const s4_acc_t *a);
double s4_acc_mae(const s4_acc_t *a);
double s4_acc_max(const s4_acc_t *a);
double s4_acc_rel(const s4_acc_t *a);

/* Rotation angle [rad] between two unit quaternions (sign-equivalent:
 * q and -q give 0). */
double s4_quat_angle(const double *qa, const double *qb);
/* Returns +1/-1 so that sgn*qb is on the same hemisphere as qa. */
int s4_quat_sign(const double *qa, const double *qb);
/* Roll/pitch/yaw [rad] (ZYX, body->NED) from scalar-first quaternion. */
void s4_quat_to_euler(const double *q, double *rpy);

#endif /* S4_METRICS_H */
