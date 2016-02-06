#include "gravity.h"

typedef struct {
    double zd;
    double r0d;
    double r0d2;
} context_t;

static void column_prepare (const context_t * context,
                    const uint16_t m, 
                    double * xk)
{
    int n = m + 1, i;
    xk[0] = (2*n - 1) * context->zd * context->r0d;
    for (i = 1, n = m + 2; n <= gravity.n; n++, i += 2) {
        xk[i] = - (n + m - 1.) / (n - m) * context->r0d2;
        xk[i+1] =  (2*n - 1.) / (n - m) * context->zd * context->r0d;
    }
}

static double potential_sum (const complex_t * cs,
                          const double * k, 
                          const uint16_t m,
                          const double * xk,
                          const complex_t * v0)
{
    double s = 0;
    complex_t v[2] = {*v0, {}};
    double sum = complex_dot_real(v[0], cs[0]);
    complex_scale_r(v[1], v[0], xk[0] * k[1], s);
    sum += complex_dot_real(v[1], cs[1]);
    int n,i,j,l;
    for (n=m+2, l=1, i=0, j=1; 
         n <= gravity.n; 
         n++, l+=2, i = (i+1) % 2, j = (j+1) % 2)
    {
        complex_add     (v[i], xk[l] * k[n-m-1], v[j], xk[l+1], s)
        complex_scale   (v[i], k[n-m], s)
        sum += complex_dot_real (v[i], cs[n-m]);
    }
    return sum;
}

static double column (const context_t * context, 
                      const uint16_t m,
                      const complex_t * cs,
                      const double * k,
                      const complex_t v0)
{
    double xk[2*(gravity.n - m - 1) + 1];
    column_prepare(context, m, xk);
    return potential_sum(cs, k, m, xk, &v0);
}

double potential (const vector_t * r) {
    complex_t tmp;
    double s = 0;
    const double rn = vector_norm(*r);
    const double r0d = gravity.r0 / rn;
    complex_t r1d = {r->x / rn, r->y / rn};
    const context_t context = {r->z / rn, r0d, r0d*r0d };
    complex_t v  = {1, 0};
    complex_scale (v, gravity.k[0], s)
    const double sum0 = column (&context, 0, &gravity.cs[0], &gravity.k[0], v);
    double sum = 0;
    int m, i;
    for (i=0, m = 1; m < gravity.n; m++) {
        i += gravity.n + 2 - m;
        // V_{n,n}
        complex_mult_r  (tmp, v, r1d)
        complex_scale_r (v, tmp, (2*m - 1) * r0d * gravity.k[i], s)
        sum += column (&context, m, &gravity.cs[i], &gravity.k[i], v);
    }
    i += 2;
    complex_mult_r  (tmp, v, r1d)
    complex_scale_r (v, tmp, (2*gravity.n - 1) * context.r0d * gravity.k[i], s)
    sum += complex_dot_real(v, gravity.cs[i]);
    return gravity.mu / rn * (sum0 + M_SQRT2 * sum);
}
