#ifndef __COSMOS_GRAVITY_H__
#define __COSMOS_GRAVITY_H__

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "complex.h"

typedef struct {
    uint16_t n;    // < 1e4
    double r0;
    double mu;
    const complex_t * cs;
    const double * k;
} gravity_t;

#define norm(x,y,z) sqrt(x*x + y*y + z*z)

double potential(const gravity_t * gravity, double x, double y, double z);

const gravity_t * gravity_get(const char * name);

#endif//__COSMOS_GRAVITY_H__
