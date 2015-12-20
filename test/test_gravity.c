#include "gravity.h"

#include <stdio.h>

static int potential44 (gravity_t * gravity, double x, double y, double z, double epsilon) {
    complex_t cs[15] = {};
    gravity->cs = cs;
    const complex_t ones = {1, 1};
    double rn = norm(x, y, z);
    double rd = gravity->r0/rn;
    double rr0[] = {1, rd, rd*rd, rd*rd*rd, rd*rd*rd*rd};
    double mu_d = gravity->mu / rn;
    double zd = z / rn;
    double zz[] = {1, zd, zd*zd, zd*zd*zd, zd*zd*zd*zd};
    complex_t r1d = {x/rn, y/rn};
    complex_t rr1[] = {{1,0}, r1d, r1d, r1d, r1d};
    complex_mult(rr1[2], rr1[1])
    complex_mult(rr1[3], rr1[2])
    complex_mult(rr1[4], rr1[3])

    double expected[] = {
        // m = 0
        mu_d * complex_dot_real(ones, rr1[0]) * rr0[0]*      zz[0],
        mu_d * complex_dot_real(ones, rr1[0]) * rr0[1]*      zz[1]                   * sqrt(3.),
        mu_d * complex_dot_real(ones, rr1[0]) * rr0[2]*(  3.*zz[2]             - 1.) * sqrt(5./4.),
        mu_d * complex_dot_real(ones, rr1[0]) * rr0[3]*( 15.*zz[3] - 9.* zz[1]     ) * sqrt(7.)/6.,
        mu_d * complex_dot_real(ones, rr1[0]) * rr0[4]*(105.*zz[4] - 90.*zz[2] + 9.) * sqrt(9.)/24.,
        // m = 1
        mu_d * complex_dot_real(ones, rr1[1]) * rr0[1]*      zz[0]              * sqrt(3.),
        mu_d * complex_dot_real(ones, rr1[1]) * rr0[2]*   3.*zz[1]              * sqrt(5./3.),
        mu_d * complex_dot_real(ones, rr1[1]) * rr0[3]*( 15.*zz[2] - 3.       ) * sqrt(7./24.),
        mu_d * complex_dot_real(ones, rr1[1]) * rr0[4]*(105.*zz[3] - 45.*zz[1]) * sqrt(9./10.)/6.,
        // m = 2
        mu_d * complex_dot_real(ones, rr1[2]) * rr0[2]*   3.*zz[0]        * sqrt(5./12.),
        mu_d * complex_dot_real(ones, rr1[2]) * rr0[3]*  15.*zz[1]        * sqrt(7./60.),
        mu_d * complex_dot_real(ones, rr1[2]) * rr0[4]*(105.*zz[2] - 15.) * sqrt(9./180.)/2.,
        // m = 3
        mu_d * complex_dot_real(ones, rr1[3]) * rr0[3]*  15.*zz[0] * sqrt(7./360.),
        mu_d * complex_dot_real(ones, rr1[3]) * rr0[4]* 105.*zz[1] * sqrt(9./2520.),
        // m = 4
        mu_d * complex_dot_real(ones, rr1[4]) * rr0[4]* 105.*zz[0] * sqrt(9./20160.)
    };
    double actual = 0;
    double delta = 0;
    int i;
    for (i=0; i<15; i++) {
        complex_set(cs[i], ones)
        actual = potential(gravity, x, y, z);
        delta  = 2 * fabs((actual - expected[i]) / (actual + expected[i]));
        if (delta > epsilon) {
            fprintf(stderr, "#%d delta = %+.12e actual = %+.12e expected = %+.12e\n", i, delta, actual, expected[i]);
            return 1;
        }
        complex_set_zero(cs[i])
    }
    return 0;
}

int test_potential() {
    int res = 0;
    double epsilon = 1e-12;

    const gravity_t * gem4 = gravity_get("GEM4");
    double k[15];
    gravity_t gravity = {4, 1., 1., NULL, k};
    int n,m, i=0, gem4_i;
    for (m=0; m<=4; m++) {
        gem4_i = (2*gem4->n + 3 - m)*m/2;
        for (n=m; n <= 4; n++) {
            fprintf(stderr, "i = %d gem4_i = %d\n", i, gem4_i);
            k[i++] = gem4->k[gem4_i++];
        }
    }
    for (i=0; i<15; i++) fprintf(stderr, "%d %f\n", i, k[i]);

    int size = 2 + 11 * 24;
    struct{double lmb,phi;} angles[size];
    double lmb, phi;
    angles[0].lmb = 0;
    angles[0].phi = - M_PI / 2;
    angles[1].lmb = 0;
    angles[1].phi = M_PI / 2;
    i = 2;
    for (lmb = 0; lmb < 2*M_PI - 0.1; lmb += M_PI/12) {
        for (phi = - 5*M_PI/12; phi < M_PI/2 - 0.1; phi += M_PI/12) {
            angles[i].lmb = lmb;
            angles[i].phi = phi;
            i++;
        }
    }
    double x,y,z;
    for (i=0; i<size; i++) {
        fprintf(stderr, " angles : %f, %f \n", angles[i].lmb, angles[i].phi);
        x = cos(angles[i].phi) * cos(angles[i].lmb);
        y = cos(angles[i].phi) * sin(angles[i].lmb);
        z = sin(angles[i].phi);
        fprintf(stderr, " r = %f, %f, %f \n", x, y, z);
        if (res = potential44(&gravity, x,y,z, epsilon)) return res;
    }
}
