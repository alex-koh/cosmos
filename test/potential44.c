#include <stdio.h>

#include "potential44.h"

static double compare (const double epsilon, double lmb, double phi) {
    double s;
    const double p0 = gravity.mu / gravity.r0;
    const double abs_epsilon = epsilon * gravity.mu / gravity.r0;
    double rn = gravity.r0*1.16; // h ~ 1000 km
    vector_t r = { cos(phi) * cos(lmb), cos(phi) * sin(lmb), sin(phi) };
    vector_scale(r, s, rn)
    double expected = potential44(&r);
    double actual = potential(&r);
    double delta = fabs((actual - expected) / (actual + expected));
    if (delta > epsilon) {
        fprintf(stderr, "ERROR epsilon = %.16e\n", epsilon);
        fprintf(stderr, "ERROR potential0 = %.16e\n", p0);
        fprintf(stderr, "ERROR abs_epsilon = %.16e\n", abs_epsilon);
        fprintf(stderr, "ERROR angles : %f, %f \n", lmb*180/M_PI, phi*180/M_PI);
        fprintf(stderr, "ERROR delta = %+.16e actual = %+.16e expected = %+.16e\n", 
                delta, actual, expected);
            return 1;
    }
    return 0;
}

#define get_type(type, fmt)                         \
    static type get_ ## type(const char * str) {    \
    type s;                                         \
    sscanf(str, fmt, &s);                           \
    return s;                                       \
}

get_type(double, "%lf")
get_type(int,    "%d")

int main(int argc, const char ** argv) {
    if (argc != 4) return 1;
    const double epsilon = get_double(argv[1]);
    const int lmb_size = get_int(argv[2]);
    const int phi_size = get_int(argv[3]);

    int res = 0;
    if (res = compare(epsilon, 0., M_PI_2)) return res;
    if (res = compare(epsilon, 0., - M_PI_2)) return res;
    int lmb, phi;
    for (lmb = 0; lmb < lmb_size; lmb++)
        for (phi = 1; phi < phi_size; phi++)
            if (res = compare(epsilon, 2*lmb*M_PI/lmb_size, 
                        - M_PI_2 + 2*phi*M_PI/phi_size)) return res;
    return 0;
}
