from sympy import Symbol, symbols, Function
from sympy import ccode, srepr, print_tree
from sympy import re, I
from sympy import Rational
from sympy import sqrt

import sys

x,y,z, mu,r0, c,s = symbols(
    "r->x r->y r->z "
    "gravity.mu gravity.r0 "
    "cs.x cs.y", real=True)

rn = sqrt(x*x+y*y+z*z)

cs = c - I * s

mu_d = mu / rn
r0d   = r0 / rn
zd   = z  / rn
r1d  = (x + I*y)/rn

n, m = int(sys.argv[1]), int(sys.argv[2])
e = {
# m = 0
     (0,0) : lambda : mu_d * re(cs),
     (1,0) : lambda : mu_d * re(cs) * r0d   *     zd                    * sqrt(3),
     (2,0) : lambda : mu_d * re(cs) * r0d**2*(  3*zd**2            - 1) * sqrt(Rational(5,4)),
     (3,0) : lambda : mu_d * re(cs) * r0d**3*( 15*zd**3 -  9*zd       ) * sqrt(7)/6,
     (4,0) : lambda : mu_d * re(cs) * r0d**4*(105*zd**4 - 90*zd**2 + 9) * sqrt(9)/24,
# m = 1
     (1,1) : lambda : mu_d * re(cs*r1d) * r0d                        * sqrt(3),
     (2,1) : lambda : mu_d * re(cs*r1d) * r0d**2*   3*zd             * sqrt(Rational(5,3)),
     (3,1) : lambda : mu_d * re(cs*r1d) * r0d**3*( 15*zd**2 -  3   ) * sqrt(Rational(7,24)),
     (4,1) : lambda : mu_d * re(cs*r1d) * r0d**4*(105*zd**3 - 45*zd) * sqrt(Rational(9,10))/6,
# m = 2
     (2,2) : lambda : mu_d * re(cs*r1d**2) * r0d**2*   3             * sqrt(Rational(5,12)),
     (3,2) : lambda : mu_d * re(cs*r1d**2) * r0d**3*  15*zd          * sqrt(Rational(7,60)),
     (4,2) : lambda : mu_d * re(cs*r1d**2) * r0d**4*(105*zd**2 - 15) * sqrt(Rational(9,180))/2,
# m = 3
     (3,3) : lambda : mu_d * re(cs*r1d**3) * r0d**3*  15    * sqrt(Rational(7,360)),
     (4,3) : lambda : mu_d * re(cs*r1d**3) * r0d**4* 105*zd * sqrt(Rational(9,2520)),
# m = 4
     (4,4) : lambda : mu_d * re(cs*r1d**4) * r0d**4* 105* sqrt(Rational(9,20160))
}[(n,m)]()

index = [[0,0,0,0,0],[1,5,0,0,0],[2,6,9,0,0],[3,7,10,12,0],[4,8,11,13,14]]

f = open(sys.argv[3], "w")
try:
    f.write("\n".join([
        "#include \"potential44.h\"",
        "",
        "double potential44(const vector_t * r) {",
        "    double " + ccode(rn, assign_to="rn"),
        "    complex_t cs = gravity.cs[%d];" % (index[n][m],),
        "    return " + ccode(e.subs(rn, Symbol("rn"))) + ";",
        "}",
        "" ,
        "vector_t gradient44(const vector_t * r) {",
        "    double " + ccode(rn, assign_to="rn"),
        "    complex_t cs = gravity.cs[%d];" % (index[n][m],),
        "    vector_t g = {",
        "       " + ccode(e.diff(x).subs(rn, Symbol("rn"))) + ",",
        "       " + ccode(e.diff(y).subs(rn, Symbol("rn"))) + ",",
        "       " + ccode(e.diff(z).subs(rn, Symbol("rn")))      ,
        "    };",
        "    return g;",
        "}"
    ]));
finally:
    f.close()
