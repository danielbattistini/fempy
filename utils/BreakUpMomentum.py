'''
Script to compute the break-up momentum of a particle into 2 daughters.
Obtained with the conservation of the quadrimomentum.


'''

import argparse
from ROOT import TDatabasePDG
import math

def ComputeSof2Particles(m1, k1, m2, k2):
    return m1*m1 + m2*m2 + 2*math.sqrt(m1*m1 + k1*k1)*math.sqrt(m2*m2 + k2*k2) + 2*k1*k2

def kaellen(m1, m2, m3):
    '''Triangular Kaellen function'''
    return m1**2 + m2**2 + m3**2 - 2*m1*m2 - 2*m2*m3 - 2*m1*m3


def BreakUpMomentum(M, m1, m2):
    '''Compute the break-up momentum p of the decay P(M, 0) -> p1(E1, p) + p2(E2, -p)'''
    return kaellen(M**2, m1**2, m2**2)**0.5/2/M

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mass', action='store_true', default=False)
    parser.add_argument('P', help='PDG code or mass of mother particle')
    parser.add_argument('p1', help='PDG code or mass of first daughter')
    parser.add_argument('p2', help='PDG code or mass of second daughter')
    args = parser.parse_args()

    if args.mass:
        M = float(args.P)
        m1 = float(args.p1)
        m2 = float(args.p2)

        print(f'Break up momentum for {M} --> {m1} {m2} ='
              f' {BreakUpMomentum(M, m1, m2)*1000:.3f} MeV/c')

    else:
        P = TDatabasePDG.Instance().GetParticle(int(args.P))
        p1 = TDatabasePDG.Instance().GetParticle(int(args.p1))
        p2 = TDatabasePDG.Instance().GetParticle(int(args.p2))

        M = P.Mass()
        m1 = p1.Mass()
        m2 = p2.Mass()

        print(f'Mass of {P.GetName()} = {M} GeV/c2')
        print(f'Mass of {p1.GetName()} = {m1} GeV/c2')
        print(f'Mass of {p2.GetName()} = {m2} GeV/c2')

        print(f'Q-value: {M - m1 - m2:.3f} GeV/c2')

        print(f'Break up momentum for {P.GetName()} --> {p1.GetName()} {p2.GetName()} = '
              f'{BreakUpMomentum(M, m1, m2)*1000:.3f} MeV/c')
