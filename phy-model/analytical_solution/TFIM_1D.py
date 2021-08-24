import numpy as np
import scipy.integrate as integrate
from numpy import cosh, cos, log, sqrt
'''

1D Transverse Field Ising Model(1D TFIM)

H = \sum_i -J*\sigma_{i}\sigma_{i+1} - h*\sigma_{i}

'''


# beta = 1/kT
# Jordan–Wigner transformation
def lnZ(beta, J, h):
    lam = J / h
    w = lambda x: (1 + 2 * lam * cos(x) + lam ** 2) ** 0.5
    a = lambda x: log(cosh(h * beta * w(x))) / np.pi
    f2 = integrate.quad(a, 0, np.pi, epsabs=1.0e-13, epsrel=1.0e-13)
    lnZ = log(2) + f2[0]
    return lnZ