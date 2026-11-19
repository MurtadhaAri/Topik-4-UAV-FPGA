/*
 * ekf_sensor_packet.h -- COPY of ICD ROLE 5, butir 1 ("sensor_packet_t").
 * Source: Kesepakatan Teknis / ICD ES-EKF.  Replace with the team's version
 * when the repository is available (same include guard).
 */
#ifndef EKF_SENSOR_PACKET_H
#define EKF_SENSOR_PACKET_H

#include <stdint.h>
#include "ekf_types.h"

typedef enum {
    SENS_IMU_ONLY = 0,     /* IMU prediction only                          */
    SENS_GNSS_POS_VEL = 1, /* GNSS position (3) + velocity (3)             */
    SENS_BARO = 2,         /* Barometer altitude (1): meas[0] = -altitude  */
    SENS_MAG = 3           /* 3-axis magnetometer (3)                      */
} sensor_type_e;

typedef struct {
    uint32_t timestamp_us;   /* time [us]                                  */
    ekf_scalar_t dt;         /* delta time since last IMU [s]              */
    sensor_type_e type;      /* which exteroceptive sensor is active       */
    ekf_scalar_t imu_acc[3]; /* [ax,ay,az] body frame [m/s^2]              */
    ekf_scalar_t imu_gyr[3]; /* [gx,gy,gz] body frame [rad/s]              */
    ekf_scalar_t meas[6];    /* exteroceptive payload according to 'type'  */
} sensor_frame_t;

#endif /* EKF_SENSOR_PACKET_H */
