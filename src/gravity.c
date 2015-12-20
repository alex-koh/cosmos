#include "gravity.h"

#include <stdio.h>

typedef struct {
    complex_t r1d;
    double zd;
    double r0d;
} context_t;

static double column (const gravity_t * gravity, 
                      const context_t * context, 
                      const uint16_t m,
                      uint32_t i,
                      complex_t v0)
{
    const double r0d2 = context->r0d * context->r0d;
//    uint32_t i = (2 * gravity->n + 3 - m) * m / 2;
    double s;
    // \bar V_{m, m}
//    complex_scale(v0, s, gravity->k[i])
    double sum = complex_dot_real(v0, gravity->cs[i]);
    // \bar V_{m + 1, m}
    int n = m + 1;
    i++;
    complex_t v1 = v0;
    complex_scale(v1, s, (2*n - 1) * context->zd * context->r0d * gravity->k[i]);
    complex_t *vp0 = &v0, *vp1 = &v1, *vpt;
    sum += complex_dot_real(v1, gravity->cs[i]);
    // \bar V_{m + k, m} (k >= 2)
    for(n = m + 2; n <= gravity->n; n++) {
        i++;
        complex_scale((*vp0), s,   - (n + m - 1.) / (n - m) * r0d2 * gravity->k[i - 1])
        complex_add((*vp0), (*vp1), s, (2*n - 1.) / (n - m) * context->zd * context->r0d)
        complex_scale((*vp0), s, gravity->k[i])
        sum += complex_dot_real((*vp0), gravity->cs[i]);
        // swap v0 and v1
        vpt = vp0; vp0 = vp1; vp1 = vpt;
    }
    return sum;
}

double potential (const gravity_t * gravity, 
                  double x, 
                  double y, 
                  double z)
{
    double s;
    const double rn = norm(x,y,z);
    const context_t context = {{x / rn, y / rn}, z / rn, gravity->r0 / rn };
    uint32_t i=0;
    complex_t v  = {1, 0};
    complex_scale (v, s, gravity->k[i])
    const double sum0 = column (gravity, &context, 0, i, v);
    double sum = 0;
    int m;
    for (m = 1; m <= gravity->n; m++) {
        i += gravity->n + 2 - m;
        // V_{n,n}
        complex_mult (v, context.r1d)
        complex_scale (v, s, (2*m - 1) * context.r0d);
        complex_scale (v, s, gravity->k[i])
        sum += column (gravity, &context, m, i, v);
    }
    fprintf(stderr, "sum0 = %f sum = %f\n", sum0, sum);
    return gravity->mu / rn * (sum0 + sqrt(2) * sum);
}
