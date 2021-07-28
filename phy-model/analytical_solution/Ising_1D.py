import TFIM_1D

'''

1D Ising model

H = \sum_i -J*\sigma_{i}\sigma_{i+1}

'''


# beta = 1/kT
# https://en.wikipedia.org/wiki/Ising_model#Ising's_exact_solution
def free_energy(J, beta):

    return TFIM_1D.free_energy(J, 0, beta)
