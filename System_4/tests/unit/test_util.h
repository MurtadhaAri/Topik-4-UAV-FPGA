/* Minimal test harness: CHECK* macros, per-binary summary, non-zero exit on failure. */
#ifndef TEST_UTIL_H
#define TEST_UTIL_H
#include <math.h>
#include <stdio.h>

static int t_checks = 0, t_fails = 0;

#define CHECK(cond, msg)                                                                      \
    do {                                                                                      \
        t_checks++;                                                                           \
        if (!(cond)) {                                                                        \
            t_fails++;                                                                        \
            printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, (msg));                           \
        }                                                                                     \
    } while (0)

#define CHECK_EQ(a, b, msg)                                                                   \
    do {                                                                                      \
        long long _a = (long long)(a), _b = (long long)(b);                                   \
        t_checks++;                                                                           \
        if (_a != _b) {                                                                       \
            t_fails++;                                                                        \
            printf("  FAIL %s:%d  %s  (got %lld, expected %lld)\n", __FILE__, __LINE__, (msg), \
                   _a, _b);                                                                   \
        }                                                                                     \
    } while (0)

#define CHECK_NEAR(a, b, tol, msg)                                                            \
    do {                                                                                      \
        double _a = (double)(a), _b = (double)(b), _t = (double)(tol);                        \
        t_checks++;                                                                           \
        if (!(fabs(_a - _b) <= _t)) {                                                         \
            t_fails++;                                                                        \
            printf("  FAIL %s:%d  %s  (got %.9g, expected %.9g, tol %.3g)\n", __FILE__,       \
                   __LINE__, (msg), _a, _b, _t);                                              \
        }                                                                                     \
    } while (0)

#define TEST_BEGIN(name) printf("== %s\n", (name))
#define TEST_END()                                                                            \
    (printf("-- %d checks, %d failed\n", t_checks, t_fails), (t_fails ? 1 : 0))

static unsigned t_rng = 0x9E3779B9u;
static inline unsigned t_rand(void)
{
    t_rng = t_rng * 1664525u + 1013904223u;
    return t_rng;
}
#endif
