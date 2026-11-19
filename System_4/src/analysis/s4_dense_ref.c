/*
 * s4_dense_ref.c -- independent dense double-precision covariance EKF used ONLY
 * to measure the error of the float32 UDU reference (analysis/test tooling).
 * Predict: P = Phi P Phi^T + diag(qc dt) ; update: sequential scalar, unit-row H.
 */
#include <math.h>
#include <string.h>

#include "ekf_s4_models.h"
#include "s4_dataset.h"
#include "s4_dense_ref.h"

/* ---------- independent dense double-precision covariance EKF (test only) ---------- */
static double P[15][15];

static void dense_predict(const double *Phi, double dt, const s4_params_t *prm)
{
    double T[15][15], N[15][15];
    int i, j, k;
    for (i = 0; i < 15; i++)
        for (j = 0; j < 15; j++) {
            double s = 0;
            for (k = 0; k < 15; k++) s += Phi[i * 15 + k] * P[k][j];
            T[i][j] = s;
        }
    for (i = 0; i < 15; i++)
        for (j = 0; j < 15; j++) {
            double s = 0;
            for (k = 0; k < 15; k++) s += T[i][k] * Phi[j * 15 + k];
            N[i][j] = s;
        }
    for (i = 0; i < 15; i++) N[i][i] += (double)prm->qc[i / 3] * dt;
    memcpy(P, N, sizeof P);
}

static void dense_update(int i0, double R)
{
    double k[15], a = P[i0][i0] + R;
    int i, j;
    for (i = 0; i < 15; i++) k[i] = P[i][i0];
    for (i = 0; i < 15; i++)
        for (j = 0; j < 15; j++) P[i][j] -= k[i] * k[j] / a;
}

static sensor_frame_t frames[S4_MAX_FRAMES];
static s4_truth_t truth[S4_MAX_FRAMES];

double s4_verify_float_vs_dense(s4_scenario_t sc, double dur)
{
    static s4f_model_t m;
    s4_params_t prm;
    int n, k, i, j, q;
    double prev = 0.0, worst = 0.0;
    float v0[3];
    s4_params_default(&prm);
    n = s4_dataset_generate(sc, dur, frames, truth, S4_MAX_FRAMES);
    s4f_init(&m, &prm);
    s4_dataset_filter_init(sc, v0);
    s4f_set_velocity(&m, v0);
    memset(P, 0, sizeof P);
    for (i = 0; i < 15; i++) P[i][i] = prm.D0[i / 3];
    for (k = 0; k < n; k++) {
        double dtc = prev + frames[k].dt;
        s4f_step(&m, &frames[k]);
        if (m.dt_acc == 0.0f && dtc > 0.0) {
            double Phi[225], a[3], w[3], R[9], Sa[9], Sw[9];
            memset(Phi, 0, sizeof Phi);
            for (i = 0; i < 15; i++) Phi[i * 15 + i] = 1.0;
            for (i = 0; i < 3; i++) { a[i] = m.a_unb[i]; w[i] = m.w_unb[i]; }
            for (i = 0; i < 9; i++) R[i] = m.R[i];
            Sa[0] = 0; Sa[1] = -a[2]; Sa[2] = a[1]; Sa[3] = a[2]; Sa[4] = 0; Sa[5] = -a[0];
            Sa[6] = -a[1]; Sa[7] = a[0]; Sa[8] = 0;
            Sw[0] = 0; Sw[1] = -w[2]; Sw[2] = w[1]; Sw[3] = w[2]; Sw[4] = 0; Sw[5] = -w[0];
            Sw[6] = -w[1]; Sw[7] = w[0]; Sw[8] = 0;
            for (i = 0; i < 3; i++) {
                Phi[i * 15 + 3 + i] = dtc;
                Phi[(6 + i) * 15 + 12 + i] = -dtc;
                for (j = 0; j < 3; j++) {
                    double ra = 0;
                    for (q = 0; q < 3; q++) ra += R[3 * i + q] * Sa[3 * q + j];
                    Phi[(3 + i) * 15 + 6 + j] = -ra * dtc;
                    Phi[(3 + i) * 15 + 9 + j] = -R[3 * i + j] * dtc;
                    Phi[(6 + i) * 15 + 6 + j] += -Sw[3 * i + j] * dtc;
                }
            }
            dense_predict(Phi, dtc, &prm);
            prev = 0.0;
        } else {
            prev = dtc;
        }
        if (frames[k].type == SENS_GNSS_POS_VEL) {
            for (i = 0; i < 3; i++) dense_update(i, prm.r_pos);
            for (i = 0; i < 3; i++) dense_update(3 + i, prm.r_vel);
        } else if (frames[k].type == SENS_BARO) {
            dense_update(2, prm.r_baro);
        }
    }
    /* reconstruct diag(U D U^T) and compare, relative */
    for (i = 0; i < 15; i++) {
        double s = 0.0;
        for (j = i; j < 15; j++) {
            double u = m.ekf.U[U_IDX(i, j)];
            s += u * u * (double)m.ekf.D[j];
        }
        if (fabs(s - P[i][i]) / P[i][i] > worst) worst = fabs(s - P[i][i]) / P[i][i];
    }
    return worst;
}

