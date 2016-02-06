#ifndef __COSMOS_COMPLEX_H__
#define __COSMOS_COMPLEX_H__

typedef struct {
    double x,y;
} complex_t;

#define complex_dot_real(z1, z2) (z1.x * z2.x + z1.y * z2.y)

#define complex_add(r, a, z, b, s) \
    s = a; \
    r.x *= s; \
    r.y *= s; \
    s = b; \
    r.x += s * z.x; \
    r.y += s * z.y; 

#define complex_scale(r, exp, s) \
    s = exp; \
    r.x *= s; \
    r.y *= s;

#define complex_scale_r(r, z, exp, s) \
    s = exp; \
    r.x = s * z.x; \
    r.y = s * z.y;

#define complex_mult_r(r, z1, z2) \
    r.x = z1.x * z2.x - z1.y * z2.y ; \
    r.y = z1.x * z2.y + z1.y * z2.x ; 

#endif//__COSMOS_COMPLEX_H__
