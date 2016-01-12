#include "gravity.h"

typedef struct {
    complex_t r1d;
    double zd;
    double r0d;
    double r0d2;
} context_t;

static double column (const gravity_t * gravity, 
                      const context_t * context, 
                      const uint16_t m,
                      uint32_t i,
                      complex_t v0)
{
    double s;
    // \bar V_{m, m}
    double sum = complex_dot_real(v0, gravity->cs[i]);
    // \bar V_{m + 1, m}
    int n = m + 1;
    i++;
    complex_t v1 = v0;
    complex_scale(v1, s, (2*n - 1) * context->zd * context->r0d * gravity->k[i]);
    complex_t *vp0 = &v0, *vp1 = &v1, *vpt;
    sum += complex_dot_real(v1, gravity->cs[i]);
    // \bar V_{m + k, m} (k >= 2)
    for (n = m + 2; n <= gravity->n; n++) {
        complex_scale((*vp0), s,   - (n + m - 1.) / (n - m) * context->r0d2 * gravity->k[i])
        complex_add((*vp0), (*vp1), s, (2*n - 1.) / (n - m) * context->zd * context->r0d)
        i++;
        complex_scale((*vp0), s, gravity->k[i])
        sum += complex_dot_real((*vp0), gravity->cs[i]);
        // swap v0 and v1
        vpt = vp0; vp0 = vp1; vp1 = vpt;
    }
    return sum;
}

double potential (const gravity_t * gravity, const vector_t * r) {
    double s;
    const double rn = vector_norm(*r);
    const double r0d = gravity->r0 / rn;
    const context_t context = {{r->x / rn, r->y / rn}, r->z / rn, r0d, r0d*r0d };
    uint32_t i=0;
    complex_t v  = {1, 0};
    complex_scale (v, s, gravity->k[i])
    const double sum0 = column (gravity, &context, 0, i, v);
    double sum = 0;
    int m;
    for (m = 1; m < gravity->n; m++) {
        i += gravity->n + 2 - m;
        // V_{n,n}
        complex_mult (v, context.r1d, s)
        complex_scale (v, s, (2*m - 1) * context.r0d * gravity->k[i])
        sum += column (gravity, &context, m, i, v);
    }
    i += 2;
    complex_mult (v, context.r1d, s)
    complex_scale (v, s, (2*gravity->n - 1) * context.r0d * gravity->k[i])
    sum += complex_dot_real(v, gravity->cs[i]);
    return gravity->mu / rn * (sum0 + M_SQRT2 * sum);
}
