#include <stdio.h>

#include "gravity.h"

static double potential44(const vector_t * r, const int i) {
    double s;
    double rn = vector_norm(*r);
    double mu_d = gravity.mu / rn;
    double rd = gravity.r0/rn;
    double zd = r->z / rn;
    complex_t r1d = {r->x/rn, r->y/rn};

    double rr0[5]; rr0[0] = 1.;
    double zz[5];  zz [0] = 1.;
    double cs_rr1[5];
    complex_t rr1 = {1., 0};
    cs_rr1[0] = complex_dot_real(gravity.cs[i], rr1);

    int j;
    for (j=1; j<5; j++) {
        zz[j] = zz[j-1]*zd;
        rr0[j] = rr0[j-1] * rd;
        complex_mult(rr1, r1d, s)
        cs_rr1[j] = complex_dot_real(rr1, gravity.cs[i]);
    }

    // m = 0
    if (i== 0) return mu_d * cs_rr1[0] * rr0[0]*      zz[0];
    if (i== 1) return mu_d * cs_rr1[0] * rr0[1]*      zz[1]                   * sqrt(3.);
    if (i== 2) return mu_d * cs_rr1[0] * rr0[2]*(  3.*zz[2]             - 1.) * sqrt(5./4.);
    if (i== 3) return mu_d * cs_rr1[0] * rr0[3]*( 15.*zz[3] - 9.* zz[1]     ) * sqrt(7.)/6.;
    if (i== 4) return mu_d * cs_rr1[0] * rr0[4]*(105.*zz[4] - 90.*zz[2] + 9.) * sqrt(9.)/24.;
    // m = 1
    if (i== 5) return mu_d * cs_rr1[1] * rr0[1]*      zz[0]              * sqrt(3.);
    if (i== 6) return mu_d * cs_rr1[1] * rr0[2]*   3.*zz[1]              * sqrt(5./3.);
    if (i== 7) return mu_d * cs_rr1[1] * rr0[3]*( 15.*zz[2] - 3.       ) * sqrt(7./24.);
    if (i== 8) return mu_d * cs_rr1[1] * rr0[4]*(105.*zz[3] - 45.*zz[1]) * sqrt(9./10.)/6.;
    // m = 2
    if (i== 9) return mu_d * cs_rr1[2] * rr0[2]*   3.*zz[0]        * sqrt(5./12.);
    if (i==10) return mu_d * cs_rr1[2] * rr0[3]*  15.*zz[1]        * sqrt(7./60.);
    if (i==11) return mu_d * cs_rr1[2] * rr0[4]*(105.*zz[2] - 15.) * sqrt(9./180.)/2.;
    // m = 3
    if (i==12) return mu_d * cs_rr1[3] * rr0[3]*  15.*zz[0] * sqrt(7./360.);
    if (i==13) return mu_d * cs_rr1[3] * rr0[4]* 105.*zz[1] * sqrt(9./2520.);
    // m = 4
    if (i==14) return mu_d * cs_rr1[4] * rr0[4]* 105.*zz[0] * sqrt(9./20160.);
    return -1;
}

static double compare (const double epsilon, double lmb, double phi, const int i) {
    double s;
    const double p0 = gravity.mu / gravity.r0;
    const double abs_epsilon = epsilon * gravity.mu / gravity.r0;
    double rn = gravity.r0*1.16; // h ~ 1000 km
    vector_t r = { cos(phi) * cos(lmb), cos(phi) * sin(lmb), sin(phi) };
    vector_scale(r, s, rn)
    double expected = potential44(&r, i);
    double actual = potential(&r);
    double delta = fabs((actual - expected) / (actual + expected));
    if (delta > epsilon) {
        fprintf(stderr, "ERROR epsilon = %.16e\n", epsilon);
        fprintf(stderr, "ERROR potential0 = %.16e\n", p0);
        fprintf(stderr, "ERROR abs_epsilon = %.16e\n", abs_epsilon);
        fprintf(stderr, "ERROR angles : %f, %f \n", lmb*180/M_PI, phi*180/M_PI);
        fprintf(stderr, "ERROR #%d delta = %+.16e actual = %+.16e expected = %+.16e\n", 
                i, delta, actual, expected);
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

find_not_zero(int i) {
    if (i>=15) return -1;
    if (fabs(gravity.cs[i].x) > 0.01) return i;
    return find_not_zero(i+1);
}

int main(int argc, const char ** argv) {
    if (argc != 4) return 1;
    const double epsilon = get_double(argv[1]);
    const int lmb_size = get_int(argv[2]);
    const int phi_size = get_int(argv[3]);
    const int i = find_not_zero(0);
    if (i < 0) return 1;

    int res = 0;
    if (res = compare(epsilon, 0., M_PI_2, i)) return res;
    if (res = compare(epsilon, 0., - M_PI_2, i)) return res;
    int lmb, phi;
    for (lmb = 0; lmb < lmb_size; lmb++)
        for (phi = 1; phi < phi_size; phi++)
            if (res = compare(epsilon, 2*lmb*M_PI/lmb_size, 
                        - M_PI_2 + 2*phi*M_PI/phi_size, i)) return res;
    return 0;
}
