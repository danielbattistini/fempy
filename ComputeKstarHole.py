import numpy as np

from ROOT import TDatabasePDG

def compute_hole(dstar_mass_thr):
    mass_dzero = TDatabasePDG.Instance().GetParticle(421).Mass()
    mass_pi = TDatabasePDG.Instance().GetParticle(211).Mass()
    mass_dstar = TDatabasePDG.Instance().GetParticle(413).Mass()

    print(mass_dstar, mass_dzero, mass_pi)

    redmass = mass_pi*mass_dzero / (mass_pi + mass_dzero)
    deltamass = dstar_mass_thr - mass_pi

    # move to D* pi ref frame
    mplus = mass_dstar**2 - (mass_dzero + mass_pi)**2
    mminus = mass_dstar**2 - (mass_dzero - mass_pi)**2
    p= (mplus * mminus)**0.5 / 2 /mass_dstar


    e_pi = mass_pi + p*p/2/mass_pi
    beta = p / (e_pi + mass_dstar)

    gamma = (1/(1-beta**2))**0.5
    return np.sqrt(2 * deltamass * redmass), p, beta, gamma




print(compute_hole(0.245))
print(compute_hole(0.18))
print(compute_hole(0.147))
print(compute_hole(0.144))
