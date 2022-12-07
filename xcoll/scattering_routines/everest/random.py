from ._everest import lib

def get_random():
    return lib.get_random()

def get_random_ruth():
    return lib.get_random_ruth()

def get_random_gauss():
    return lib.get_random_gauss()

def set_random_seed(seed=0):
    lib.set_random_seed(seed)

def get_random_seed():
    return lib.get_random_seed()

def get_n_sampled():
    return lib.get_n_sampled()

def set_rutherford_iterations(n_iter):
    lib.set_rutherford_iterations(n_iter)

def get_rutherford_iterations():
    return lib.get_rutherford_iterations()

def set_rutherford_parameters(zatom, emr, hcut):
    lib.set_rutherford_parameters(zatom, emr, hcut)


