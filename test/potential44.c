#include <stdio.h>

#include "gravity.h"

static int potential44 (gravity_t * gravity, double lmb, double phi, double epsilon) {
    double rn = gravity->r0*1.16; // h ~ 1000 km
    double s = 0;
    vector_t r = { cos(phi) * cos(lmb), cos(phi) * sin(lmb), sin(phi) };
    vector_scale(r, s, rn)
    complex_t cs[15] = {};
    gravity->cs = cs;
    const complex_t ones = {1, 1};
    double rd = gravity->r0/rn;
    double rr0[] = {1, rd, rd*rd, rd*rd*rd, rd*rd*rd*rd};
    double mu_d = gravity->mu / rn;
    double zd = r.z / rn;
    double zz[] = {1, zd, zd*zd, zd*zd*zd, zd*zd*zd*zd};
    complex_t r1d = {r.x/rn, r.y/rn};
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
    /*
     * Because I use {@code ones} instead actual {@code cs} it is a 
     * catastrophic cancellation appeare in case when \lambda == \pi/4 + \pi*n,
     * but because results are normilized, I just ignore this delta if it < \epsilon.
     * 
     * */
    epsilon *= gravity->mu/gravity->r0;
    double actual = 0;
    double delta = 0;
    int i;
    for (i=0; i<15; i++) {
        complex_set(cs[i], ones)
        actual = potential(gravity, &r);
        if ((expected[i] > epsilon) &&  
            ((delta  = 2 * fabs((actual - expected[i]) / (actual + expected[i]))) > epsilon))
        {
            fprintf(stderr, " angles : %f, %f \n", lmb*180/M_PI, phi*180/M_PI);
            fprintf(stderr, "#%d delta = %+.16e actual = %+.16e expected = %+.16e\n", i, delta, actual, expected[i]);
            return 1;
        }
        complex_set_zero(cs[i])
    }
    return 0;
}

int main(int argc, const char ** argv) {
    fprintf(stderr, "START TESTS : potential44\n");

    double epsilon = 1e-14;
    int res = 0;

    gravity_t test_gravity = gravity;

    int size = 2 + 11 * 24;
    struct{double lmb,phi;} angles[size];
    double lmb, phi;
    angles[0].lmb = 0;
    angles[0].phi = - M_PI / 2;
    angles[1].lmb = 0;
    angles[1].phi = M_PI / 2;
    int i = 2;
    for (lmb = 0; lmb < 2*M_PI - 0.1; lmb += M_PI/12) {
        for (phi = - 5*M_PI/12; phi < M_PI/2 - 0.1; phi += M_PI/12) {
            angles[i].lmb = lmb;
            angles[i].phi = phi;
            i++;
        }
    }
    // normal Earth's constants
    test_gravity.r0 = 6.4e6;
    test_gravity.mu = 4.0e14;
    for (i=0; i<size; i++) if (res = potential44(&test_gravity, angles[i].lmb, angles[i].phi, epsilon)) return res;

    // abnormal constants
    test_gravity.r0 = 6.4;
    test_gravity.mu = 4.0e25;
    for (i=0; i<size; i++) if (res = potential44(&test_gravity, angles[i].lmb, angles[i].phi, epsilon)) return res;

    fprintf(stderr, "SUCCESS : potential44");
    return 0;
}
