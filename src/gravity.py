# -*- coding: utf-8 -*-

"""
    Calculate gravity potential in form of Холщевников 2005 algorithm

    todo rewrite to OpenCL
    todo move this realization to tests
    todo create function to calc n-th derevative by r of potential

"""


import numpy as np


# conjugat cs

class Gravity(object):
    def __init__(self, n, cs, mu, r0):
        self.n = n;
        if n < 1: raise Exception("wrong n")
        size = (n+2)*(n+1)/2
        self.size = size
        # create index to last elements in line [0, 2, 5, 9, 14, ...]
        last_idx = np.arange(1, n+2)
        last_idx[0] = 0
        last_idx.cumsum(out=last_idx)
        self.last_idx = last_idx
        # index to penultimate elements in line [1, 4, 8, 13, ...]
        last_idx = last_idx[1:]
        self.last_m1_idx = last_idx - np.ones_like(last_idx)
        # line indexes [[0], [1,2], [3,4,5], [6,7,8,9], ...]
        self.idx=[]
        self.r0 = r0
        b = 0
        for e in np.arange(1,n+2).cumsum():
            self.idx.append(range(b,e))
            b = e
        self.cs = cs.conj()[:size]
        self.mu = mu
        # todo this formuls is not full normalized
        K = np.ndarray((size,), dtype=np.float64)
        K[0] = 1
        for n in range(1,self.n + 1):
            idx = self.idx[n]
            K[idx[0]] = 1/np.sqrt(2*n+1)
            K[idx[1]] = np.sqrt(n*(n+1.)/2./(2*n+1))
            for m in range(2, n + 1):
                K[idx[m]] = K[idx[m-1]] * np.sqrt(n**2 - m**2 + n + m)
        self.K = K
        # todo =====

    def calc_ab(self, cs):
        return cs


    def calc_v(self, r):
        if r is  None: raise Exception("r is none")
        # inner conteiner of v-coeffs
        v = np.ndarray((self.size, ), dtype=np.complex128)
        x,y,z = r
        # calc degrees of module of r-vector r
        rn = np.linalg.norm(r)
        r2 = rn*rn
        
        # calc last coeffs in all lines
        v_last = np.arange(-1, 2 * self.n + 1, 2, dtype=np.complex128)
        v_last[0] = 1./rn
        v_last[1:] *= (x+1j*y) / r2 * self.r0 # todo r0^n i move r0 here to exclude overflow
        v_last.cumprod(out=v_last)
        v[self.last_idx] = v_last
        
        # calc last - 1 coefs
        v_m1_last = np.arange(1, 2 * self.n + 1, 2, dtype=np.complex128)
        v_m1_last *= z / r2 * self.r0 # todo r0^n
        v_m1_last *= v_last[:-1]
        v[self.last_m1_idx] = v_m1_last
        
        # calc main coefs for all lines
        for n in range(2, self.n+1):
            t = (2*n - 1) * z * v[self.idx[n-1][:-1]] * self.r0  # todo r0^n
            t -=  np.arange(n-1, 2*n-2) * v[self.idx[n-2]] * self.r0**2 # todo r0^n
            t *= 1. / r2 / (n - np.arange(n-1))
            v[self.idx[n][:-2]] = t
        v /= self.K
        return v
    
    def potential(self, r):
        v = self.calc_v(r)
        return self.mu * np.real((self.calc_v(r) * self.cs).sum())

