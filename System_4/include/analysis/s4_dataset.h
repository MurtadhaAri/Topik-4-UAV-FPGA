/*
 * s4_dataset.h -- deterministic SYNTHETIC datasets (no real flight data).
 *
 * Output uses the ICD sensor_frame_t.  Frame k carries the IMU sample for the
 * interval ((k-1)dt, k*dt]; every 100th frame also carries GNSS pos+vel
 * (5 Hz), else every 25th frame a barometer sample (20 Hz).  Rates and all
 * IMU noise levels / true biases are PROPOSED; GNSS/baro sigma follow the ICD R.
 */
#ifndef S4_DATASET_H
#define S4_DATASET_H

#include "ekf_sensor_packet.h"

typedef enum {
    S4_SC_STATIC = 0,     /* Test 1: ideal IMU + GNSS at origin (ICD sanity case) */
    S4_SC_CONST_VEL = 1,  /* Test 2 */
    S4_SC_DYNAMIC = 2,    /* Test 3 */
    S4_SC_HIGH_RANGE = 3, /* Test 4: initial errors at ~95% of the ICD bounds */
    S4_SC_COUNT
} s4_scenario_t;

#define S4_IMU_HZ 500
#define S4_IMU_DT 0.002
#define S4_MAX_FRAMES 30001

typedef struct {
    double p[3], v[3], q[4]; /* truth at the end of the frame interval */
} s4_truth_t;

const char *s4_scenario_name(s4_scenario_t sc);

/* Initial nominal velocity given to BOTH models for the scenario (GNSS-style
 * warm start for constvel/dynamic; zero = ICD baseline cold start otherwise,
 * deliberately wrong for highrange). */
void s4_dataset_filter_init(s4_scenario_t sc, float v0[3]);

/* Fills frames[0..n-1] and truth[0..n-1]; returns n (<= max_frames). */
int s4_dataset_generate(s4_scenario_t sc, double duration_s, sensor_frame_t *frames,
                        s4_truth_t *truth, int max_frames);

#endif /* S4_DATASET_H */
