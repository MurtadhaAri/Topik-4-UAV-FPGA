/*
 * ekf_types.h -- COPY of ICD section 5 ("Kontrak Header Bahasa C").
 *
 * Source : Kesepakatan_Teknis_ICD_ES-EKF.docx, Bagian I, butir 5.
 * Status : transcribed by System 4 because no team repository checkout was
 *          available.  When the team repository is available, put ITS
 *          ekf_types.h earlier on the include path; this file is guarded by
 *          the same EKF_TYPES_H macro and will then be skipped.
 * Do NOT edit the contract here.
 */
#ifndef EKF_TYPES_H
#define EKF_TYPES_H

#include <stdint.h>
#include <stdbool.h>

#define NOM_DIM 16
#define ERR_DIM 15
#define U_UPPER_SIZE 120 /* ERR_DIM * (ERR_DIM + 1) / 2 */

/* Index of the error state (dx[15]) and row/column of the 15x15 matrices */
#define IDX_DP_X 0
#define IDX_DP_Y 1
#define IDX_DP_Z 2
#define IDX_DV_X 3
#define IDX_DV_Y 4
#define IDX_DV_Z 5
#define IDX_DTH_X 6
#define IDX_DTH_Y 7
#define IDX_DTH_Z 8
#define IDX_DBA_X 9
#define IDX_DBA_Y 10
#define IDX_DBA_Z 11
#define IDX_DBG_X 12
#define IDX_DBG_Y 13
#define IDX_DBG_Z 14

/* Standard row-major 15x15 addressing */
#define MAT15_IDX(r, c) ((r) * ERR_DIM + (c))

/* Upper-triangular column-major addressing for U of UDU^T. Requires r <= c. */
#define U_IDX(r, c) (((c) * ((c) + 1)) / 2 + (r))

typedef float ekf_scalar_t; /* Orang 4 changes this while testing fixed-point */

typedef struct {
    ekf_scalar_t p[3];  /* Position NED [m]        */
    ekf_scalar_t v[3];  /* Velocity NED [m/s]      */
    ekf_scalar_t q[4];  /* Quaternion [qw,qx,qy,qz]*/
    ekf_scalar_t ba[3]; /* Accel bias [m/s^2]      */
    ekf_scalar_t bg[3]; /* Gyro bias [rad/s]       */
} nominal_state_t;

typedef struct {
    nominal_state_t x_nom;
    ekf_scalar_t dx[ERR_DIM];     /* Error state vector (15x1)                 */
    ekf_scalar_t U[U_UPPER_SIZE]; /* Unit upper-triangular (main diagonal 1.0) */
    ekf_scalar_t D[ERR_DIM];      /* Diagonal of D (15x1)                      */
} es_ekf_t;

#endif /* EKF_TYPES_H */
