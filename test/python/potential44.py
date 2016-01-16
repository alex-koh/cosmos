from sympy import Symbol, symbols, Function
from sympy import ccode, srepr, print_tree
from sympy import re, I
from sympy import Rational
from sympy import sqrt

x,y,z, mu,r0, c,s = symbols(
    "r.x r.y r.z "
    "gravity.mu gravity.r0 "
    "cs.x cs.y", real=True)

rn = sqrt(x*x+y*y+z*z)

cs = c - I * s

mu_d = mu / rn
r0d   = r0 / rn
zd   = z  / rn
r1d  = (x + I*y)/rn

expected = [
# m = 0
     mu_d * re(cs),
     mu_d * re(cs) * r0d   *     zd                    * sqrt(3),
     mu_d * re(cs) * r0d**2*(  3*zd**2            - 1) * sqrt(Rational(5,4)),
     mu_d * re(cs) * r0d**3*( 15*zd**3 -  9*zd       ) * sqrt(7)/6,
     mu_d * re(cs) * r0d**4*(105*zd**4 - 90*zd**2 + 9) * sqrt(9)/24,
# m = 1
     mu_d * re(cs*r1d) * r0d                        * sqrt(3),
     mu_d * re(cs*r1d) * r0d**2*   3*zd             * sqrt(Rational(5,3)),
     mu_d * re(cs*r1d) * r0d**3*( 15*zd**2 -  3   ) * sqrt(Rational(7,24)),
     mu_d * re(cs*r1d) * r0d**4*(105*zd**3 - 45*zd) * sqrt(Rational(9,10))/6,
# m = 2
     mu_d * re(cs*r1d**2) * r0d**2*   3             * sqrt(Rational(5,12)),
     mu_d * re(cs*r1d**2) * r0d**3*  15*zd          * sqrt(Rational(7,60)),
     mu_d * re(cs*r1d**2) * r0d**4*(105*zd**2 - 15) * sqrt(Rational(9,180))/2,
# m = 3
     mu_d * re(cs*r1d**3) * r0d**3*  15    * sqrt(Rational(7,360)),
     mu_d * re(cs*r1d**3) * r0d**4* 105*zd * sqrt(Rational(9,2520)),
# m = 4
     mu_d * re(cs*r1d**4) * r0d**4* 105* sqrt(Rational(9,20160))
]

for e in expected: print ccode(e.subs(rn, Symbol("rn"))), "\n"
