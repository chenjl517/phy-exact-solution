import numpy as np
import scipy.integrate as integrate
from scipy.misc import derivative


# https://en.wikipedia.org/wiki/Ising_model#Ising's_exact_solution
def free_energy(beta):
    k = 1 / np.sinh(2 * beta) ** 2

    a = lambda theta: np.log(np.cosh(2 * beta) ** 2 + np.sqrt(1 + k ** 2 - 2 * k * np.cos(theta)) / k)

    intResult, err = integrate.quad(a, 0, np.pi, epsabs=1.0e-13, epsrel=1.0e-13)

    lnZ_Ns = np.log(2) / 2 + intResult / (2 * np.pi)

    return -lnZ_Ns / beta


# U = \frac{\partial}{\partial \beta} (F \beta)
def internal_energy(beta):
    F_beta = lambda beta: free_energy(beta) * beta

    U = derivative(F_beta, beta, dx=1e-7)

    return U


def magnetization(beta):
    beta_c = np.log(1 + np.sqrt(2)) / 2

    if (beta < beta_c):
        return 0
    else:
        return (1 - np.sinh(2 * beta) ** (-4)) ** (1 / 8)
