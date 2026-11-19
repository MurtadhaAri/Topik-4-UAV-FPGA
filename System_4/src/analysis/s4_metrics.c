#include "s4_metrics.h"

#include <math.h>

void s4_acc_init(s4_acc_t *a)
{
    a->n = 0;
    a->sum_sq = a->sum_abs = a->max_abs = a->sum_rel = 0.0;
    a->n_rel = 0;
}

void s4_acc_add(s4_acc_t *a, double ref, double test)
{
    double e = fabs(test - ref);
    a->n++;
    a->sum_sq += e * e;
    a->sum_abs += e;
    if (e > a->max_abs) a->max_abs = e;
    if (fabs(ref) > S4_REL_EPS) {
        a->sum_rel += e / fabs(ref);
        a->n_rel++;
    }
}

double s4_acc_rmse(const s4_acc_t *a) { return a->n ? sqrt(a->sum_sq / (double)a->n) : 0.0; }
double s4_acc_mae(const s4_acc_t *a) { return a->n ? a->sum_abs / (double)a->n : 0.0; }
double s4_acc_max(const s4_acc_t *a) { return a->max_abs; }
double s4_acc_rel(const s4_acc_t *a) { return a->n_rel ? a->sum_rel / (double)a->n_rel : 0.0; }

int s4_quat_sign(const double *qa, const double *qb)
{
    double d = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3];
    return (d < 0.0) ? -1 : 1;
}

double s4_quat_angle(const double *qa, const double *qb)
{
    /* r = conj(qa) * qb ; angle = 2*atan2(|r_v|, |r_w|) */
    double w = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3];
    double x = qa[0] * qb[1] - qa[1] * qb[0] - qa[2] * qb[3] + qa[3] * qb[2];
    double y = qa[0] * qb[2] + qa[1] * qb[3] - qa[2] * qb[0] - qa[3] * qb[1];
    double z = qa[0] * qb[3] - qa[1] * qb[2] + qa[2] * qb[1] - qa[3] * qb[0];
    return 2.0 * atan2(sqrt(x * x + y * y + z * z), fabs(w));
}

void s4_quat_to_euler(const double *q, double *rpy)
{
    double w = q[0], x = q[1], y = q[2], z = q[3];
    double sp = 2.0 * (w * y - z * x);
    if (sp > 1.0) sp = 1.0;
    if (sp < -1.0) sp = -1.0;
    rpy[0] = atan2(2.0 * (w * x + y * z), 1.0 - 2.0 * (x * x + y * y));
    rpy[1] = asin(sp);
    rpy[2] = atan2(2.0 * (w * z + x * y), 1.0 - 2.0 * (y * y + z * z));
}
