#ifndef __COSMOS_GRAVITY_H__
#define __COSMOS_GRAVITY_H__

#include <stdint.h>

#include "complex.h"
#include "vector.h"

typedef struct {
    uint16_t n;    // < 1e4
    double r0;
    double mu;
    const complex_t * cs;
    const double * k;
} gravity_t;

extern const gravity_t gravity;

double potential(const vector_t * r);

#endif//__COSMOS_GRAVITY_H__
