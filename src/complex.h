#ifndef __COSMOS_COMPLEX_H__
#define __COSMOS_COMPLEX_H__

typedef struct {
    double x,y;
} complex_t;

#define complex_dot_real(z1, z2) (z1.x * z2.x + z1.y * z2.y)

#define complex_set_zero(z1) \
    z1.x = 0.; \
    z1.y = 0.;

#define complex_set(z1, z2) \
    z1.x = z2.x; \
    z1.y = z2.y;

#define complex_add(z1, z2, t, exp) \
    t = exp; \
    z1.x += t * z2.x; \
    z1.y += t * z2.y; 

#define complex_scale(z, t, exp) \
    t = exp; \
    z.x *= t; \
    z.y *= t;

#define complex_mult(z1, z2, t) \
    t = z1.x * z2.x - z1.y * z2.y ; \
    z1.y = z1.x * z2.y + z1.y * z2.x ; \
    z1.x = s;

#endif//__COSMOS_COMPLEX_H__
