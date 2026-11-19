/* Dense double-precision cross-check of the float32 UDU kernels (analysis tooling). */
#ifndef S4_DENSE_REF_H
#define S4_DENSE_REF_H

#include "s4_dataset.h"

/* Runs Model A and an independent dense double covariance EKF on the same
 * frames; returns the worst relative error of diag(U D U^T) at the end. */
double s4_verify_float_vs_dense(s4_scenario_t sc, double duration_s);

#endif
