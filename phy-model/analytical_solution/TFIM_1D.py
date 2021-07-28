import numpy as np

'''

1D Transverse Field Ising Model(1D TFIM)

H = \sum_i -J*\sigma_{i}\sigma_{i+1} - h*\sigma_{i}

'''


# beta = 1/kT
# https://en.wikipedia.org/wiki/Ising_model#Ising's_exact_solution
def free_energy(J, h, beta):
    free_energy = -1 / beta * np.log(np.exp(beta * J) * np.cosh(beta * h) +
                                     np.sqrt(np.exp(2 * beta * J) * np.sinh(beta * h) ** 2 + np.exp(-2 * beta * J)))
    return free_energy
