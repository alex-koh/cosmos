#ifndef __COSMOS_VECTOR_H__
#define __COSMOS_VECTOR_H__

#include <math.h>

typedef struct { double x,y,z; } vector_t;

#define vector_norm(r) sqrt((r).x*(r).x + (r).y*(r).y + (r).z*(r).z)

#define vector_scale(r, s, expr) s = expr; (r).x *= s; (r).y *= s; (r).z *= s;

#endif//__COSMOS_GRAVITY_H__
