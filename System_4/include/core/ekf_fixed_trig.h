/*
 * ekf_fixed_trig.h -- fixed-point CORDIC (rotation mode, shift-and-add).
 *
 * Angle input  : Q3.29 radians, int32_t.  Valid range [-pi, pi] (and, because
 *                the format itself spans [-4,4), values in (pi,4) and
 *                [-4,-pi) are wrapped by +-2*pi once).  Convert with
 *                cordic_angle_from_double() / fxp_from_double_q(x,32,29).
 * Output       : Q2.30, sin and cos in [-1, +1]  (+1.0 = 1<<30 exactly
 *                representable, no wrap-around).
 * Algorithm    : quadrant reduction to [-pi/2, pi/2] using
 *                sin(pi-a)=sin(a), cos(pi-a)=-cos(a); then N micro-rotations
 *                x' = x - d*(y>>i), y' = y + d*(x>>i), z' = z - d*atan(2^-i)
 *                with d = sign(z).  Initial x = K (= prod 1/sqrt(1+2^-2i),
 *                0.607252935...) pre-compensates the CORDIC gain, so no
 *                final multiplication is needed.
 * Iterations   : fixed (CORDIC_ITERATIONS, default 30, max 30).
 */
#ifndef EKF_FIXED_TRIG_H
#define EKF_FIXED_TRIG_H

#include <stdint.h>
#include "ekf_fixed_config.h"

void cordic_sincos(int32_t angle, int32_t *sin_out, int32_t *cos_out);
void cordic_sincos_n(int32_t angle, uint8_t iterations, int32_t *sin_out, int32_t *cos_out);

/* Constants exposed for tests / documentation */
extern const int32_t CORDIC_ATAN_Q29[32];
#define CORDIC_K_Q30 652032874 /* prod_{i<30} 1/sqrt(1+2^-2i) in Q2.30 */
#define CORDIC_PI_Q29 1686629713
#define CORDIC_HALFPI_Q29 843314857
#define CORDIC_TWOPI_Q29 3373259426LL

int32_t cordic_angle_from_double(double rad); /* I/O helper, uses double */

#endif /* EKF_FIXED_TRIG_H */
