#include <stdio.h>

#include "gravity.h"

static int potential44 (gravity_t * gravity, double lmb, double phi, double epsilon) {
    double x = cos(phi) * cos(lmb), y = cos(phi) * sin(lmb), z = sin(phi);
    complex_t cs[15] = {};
    gravity->cs = cs;
    const complex_t ones = {1, 1};
    double s = 0;
    double rn = norm(x, y, z);
    double rd = gravity->r0/rn;
    double rr0[] = {1, rd, rd*rd, rd*rd*rd, rd*rd*rd*rd};
    double mu_d = gravity->mu / rn;
    double zd = z / rn;
    double zz[] = {1, zd, zd*zd, zd*zd*zd, zd*zd*zd*zd};
    complex_t r1d = {x/rn, y/rn};
    complex_t rr1[] = {{1,0}, r1d, r1d, r1d, r1d};
    complex_mult(rr1[2], rr1[1], s)
    complex_mult(rr1[3], rr1[2], s)
    complex_mult(rr1[4], rr1[3], s)

    double expected[] = {
    // m = 0
        /* 0*/ mu_d * complex_dot_real(ones, rr1[0]) * rr0[0]*      zz[0],
        /* 1*/ mu_d * complex_dot_real(ones, rr1[0]) * rr0[1]*      zz[1]                   * sqrt(3.),
        /* 2*/ mu_d * complex_dot_real(ones, rr1[0]) * rr0[2]*(  3.*zz[2]             - 1.) * sqrt(5./4.),
        /* 3*/ mu_d * complex_dot_real(ones, rr1[0]) * rr0[3]*( 15.*zz[3] - 9.* zz[1]     ) * sqrt(7.)/6.,
        /* 4*/ mu_d * complex_dot_real(ones, rr1[0]) * rr0[4]*(105.*zz[4] - 90.*zz[2] + 9.) * sqrt(9.)/24.,
    // m = 1
        /* 5*/ mu_d * complex_dot_real(ones, rr1[1]) * rr0[1]*      zz[0]              * sqrt(3.),
        /* 6*/ mu_d * complex_dot_real(ones, rr1[1]) * rr0[2]*   3.*zz[1]              * sqrt(5./3.),
        /* 7*/ mu_d * complex_dot_real(ones, rr1[1]) * rr0[3]*( 15.*zz[2] - 3.       ) * sqrt(7./24.),
        /* 8*/ mu_d * complex_dot_real(ones, rr1[1]) * rr0[4]*(105.*zz[3] - 45.*zz[1]) * sqrt(9./10.)/6.,
    // m = 2
        /* 9*/ mu_d * complex_dot_real(ones, rr1[2]) * rr0[2]*   3.*zz[0]        * sqrt(5./12.),
        /*10*/ mu_d * complex_dot_real(ones, rr1[2]) * rr0[3]*  15.*zz[1]        * sqrt(7./60.),
        /*11*/ mu_d * complex_dot_real(ones, rr1[2]) * rr0[4]*(105.*zz[2] - 15.) * sqrt(9./180.)/2.,
    // m = 3
        /*12*/ mu_d * complex_dot_real(ones, rr1[3]) * rr0[3]*  15.*zz[0] * sqrt(7./360.),
        /*13*/ mu_d * complex_dot_real(ones, rr1[3]) * rr0[4]* 105.*zz[1] * sqrt(9./2520.),
    // m = 4
        /*14*/ mu_d * complex_dot_real(ones, rr1[4]) * rr0[4]* 105.*zz[0] * sqrt(9./20160.)
    };
    double actual = 0;
    double delta = 0;
    int i;
    for (i=0; i<15; i++) {
        complex_set(cs[i], ones)
        actual = potential(gravity, x, y, z);
        if ((expected[i] > epsilon) &&  ((delta  = 2 * fabs((actual - expected[i]) / (actual + expected[i]))) > epsilon)) {
            fprintf(stderr, " angles : %f, %f \n", lmb, phi);
            fprintf(stderr, "#%d delta = %+.12e actual = %+.12e expected = %+.12e\n", i, delta, actual, expected[i]);
            return 1;
        }
        complex_set_zero(cs[i])
    }
    return 0;
}

int main(int argc, const char ** argv) {
    fprintf(stderr, "START TESTS : potential44\n");

    int res = 0;
    double epsilon = 1e-12;

    double k[15];
    gravity_t test_gravity = {4, 1., 1., NULL, k};
    int n, m, i=0, gravity_i;
    for (m=0; m<=4; m++) {
        gravity_i = (2*gravity.n + 3 - m)*m/2;
        for (n=m; n <= 4; n++) {
            fprintf(stderr, "i = %d gravity_i = %d\n", i, gravity_i);
            k[i++] = gravity.k[gravity_i++];
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
    for (i=0; i<size; i++) if (res = potential44(&test_gravity, angles[i].lmb, angles[i].phi, epsilon)) return res;

    fprintf(stderr, "SUCCESS : potential44");
    return 0;
}
