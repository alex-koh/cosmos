import unittest

import numpy as np
import sys

from gdf_reader import *
from gfc_reader import *

from gravity import Gravity

def potential(cs, r0, mu, r):
    rn = np.linalg.norm(r)
    x,y,z = r
    r1 = x+1j*y

    v =  np.array([
        # n = 0
        cs[ 0]      *1./rn ,
        # n = 1
        cs[ 1]*r0   *z   /rn**3   * np.sqrt(3.),
        cs[ 2]*r0   *  r1/rn**3   * np.sqrt(3.),
        # n = 2
        cs[ 3]*r0**2*(3.*z**2       - rn**2) / (2.*rn**5)   * np.sqrt(5.),
        cs[ 4]*r0**2* 3.*z   *r1             /     rn**5    * np.sqrt(5./3.),
        cs[ 5]*r0**2* 3.     *r1**2          /     rn**5    * np.sqrt(5./12.),
        # n = 3
        cs[ 6]*r0**3*(15.*z**2       -  9.*rn**2) * z / (6.*rn**7)   * np.sqrt(7.),
        cs[ 7]*r0**3*(15.*z**2       -  3.*rn**2) *r1 / (2.*rn**7)   * np.sqrt(7./6.),
        cs[ 8]*r0**3* 15.*z   *r1**2                  /     rn**7    * np.sqrt(7./60.),
        cs[ 9]*r0**3* 15.     *r1**3                  /     rn**7    * np.sqrt(7./360.),
        # n = 4
        cs[10]*r0**4*(105.*z**4 - 90.*rn**2*z**2 + 9.*rn**4) / (24.*rn**9)   * np.sqrt(9.   ),
        cs[11]*r0**4*(105.*z**2 - 45.*rn**2) * z * r1 / (6.*rn**9)           * np.sqrt(9./10.),
        cs[12]*r0**4*(105.*z**2 - 15.*rn**2) *  r1**2 / (2.*rn**9)           * np.sqrt(9./180.),
        cs[13]*r0**4* 105.*z                 *  r1**3 /     rn**9            * np.sqrt(9./2520.),
        cs[14]*r0**4* 105.                   *  r1**4 /     rn**9            * np.sqrt(9./20160.)
    ], dtype=np.complex128)
    return mu*v.real.sum()

def geoid_b(a, f):
    b = a*(1-f)
    def func(lmb, phi, h):
        return np.array([
            a**2 * np.cos(lmb) * np.cos(phi),
            a**2 * np.sin(lmb) * np.cos(phi),
            b**2 * np.sin(phi)
        ], dtype=np.float64) / np.sqrt((a*np.cos(phi))**2 + (b*np.sin(phi))**2)
    return func

class GravityVTest(unittest.TestCase):
    def test_v(self):
        r0 = 7000000.
        cs = np.ones(15)
        g = Gravity(4, cs, 1., 1)

        size = 1000
        actual = np.ndarray((size,), dtype=np.float64)
        expected = np.ndarray((size,), dtype=np.float64)

        for i in range(size):
            r = np.array([
                r0*(1 + 10*np.random.rand()),
                r0*(1 + 10*np.random.rand()),
                r0*(1 + 10*np.random.rand())
            ], dtype=np.float64)
            actual[i] = g.potential(r)
            expected[i] = potential(cs, 1, 1, r)

        delta = np.abs(2 * (expected - actual)/(expected + actual))
        flags = delta < 1e-14

        if not flags.all():
            self.fail(msg="\n".join("%2d\t%+.10e\t%+.10e\t%+.10e" % 
                    (i, actual[i], expected[i], delta[i])
                    for i in np.arange(len(delta))[~flags]
                )
            )

    def test_gem4_potential(self):
        gdf = read_gdf(os.path.join(sys.__resources__, "gem4", "gem4-9280.gdf"))
        gfc = read_gfc(os.path.join(sys.__resources__, "gem4", "gem4.gfc"))

        g = Gravity(
            gfc["max_degree"], 
            gfc["data"], 
            gfc["earth_gravity_constant"], 
            gfc["radius"]
        )

        el = geoid_b(gdf["radiusrefpot"], gdf["flatrefpot"])

        DEG = np.pi/180.
        actual = np.ndarray(gdf["data"].shape[0], dtype=np.float64)
        h = gdf["height_over_ell"]
        i = 0
        for p in gdf["data"]:
            actual[i] = g.potential(el(p[0]*DEG, p[1]*DEG, h))
            i += 1

        delta = np.abs(2*(actual-gdf["data"][:,2])/(actual+gdf["data"][:,2]))
        flags = delta < 2e-14 #todo  2
        if not flags.all():
            self.fail(msg="\n".join("%2d\t%+.10e\t%+.10e\t%+.10e" % 
                    (i, actual[i], expected[i], delta[i])
                    for i in np.arange(len(delta))[~flags]
                )
            )

    def _test_egm96_potential(self):
        self.fail("todo. It is very very slow. I am going to use OpenCL to speed up my program")

        gdf = read_gdf(os.path.join(sys.__resources__, "egm96", "egm96-9287.gdf"))
        gfc = read_gfc(os.path.join(sys.__resources__, "egm96", "egm96.gfc"))

