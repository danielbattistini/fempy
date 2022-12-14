'''
Utility functions and helpers.
'''

import sys

from ROOT import TDatabasePDG


class Pair:
    def __init__(self, pair):
        if pair == 'DstarPi':
            self.name = 'DstarPi'
            self.mass_regions = ['sgn', 'sbr']
            self.hpdg = 413
            self.lpdg = 211
            self.max_kstar = 3000  # MeV/c
            self.norm_range = [1500, 2000]

        elif pair == 'DPi':
            self.name = 'DPi'
            self.mass_regions = ['sgn', 'sbl',  'sbr']
            self.hpdg = 411
            self.lpdg = 211
            self.max_kstar = 3000  # MeV/c
            self.norm_range = [1000, 1500]

        elif pair == 'DK':
            self.name = 'DK'
            self.mass_regions = ['sgn', 'sbl',  'sbr']
            self.hpdg = 411
            self.lpdg = 321
            self.max_kstar = 3000  # MeV/c
            self.norm_range = [1500, 2000]
        else:
            print("Error: pair not implemented. Exit!")
            sys.exit()


def is_mass_selected(mass, pt, pdg=413, selection='sgn', nsigma_mass=2., nsigma_offset=5., sideband_width=0.2, lower_Dstar_removal=1.992, upper_Dstar_removal=2.028) -> bool:
    '''
    function to perform the pT dependent selection of D+, D* candidates in
    different regions of the invariant mass distribution. Based on
    AliAnalysisTaskCharmingFemto::MassSelection.
    '''

    if selection == 'any':
        return True

    # mass shift observed in all Run2 data samples for all D-meson species
    mass_mean = TDatabasePDG.Instance().GetParticle(pdg).Mass() + 0.0025
    if pdg == 411:
        massWidth = 0.006758 + pt * 0.0005124
    elif pdg == 413:
        mDstarPDG = TDatabasePDG.Instance().GetParticle(413).Mass()
        mD0PDG = TDatabasePDG.Instance().GetParticle(421).Mass()
        mass_mean = mDstarPDG-mD0PDG  # no extra mass shift because it is deltamass
        massWidth = 0.00124673 - pt * 0.000340426 + pt * pt * 4.40729e-05
        if pt > 4 and pt < 5:
            massWidth = 0.00104329 - 0.000113275 * pt
        elif pt >= 5:
            massWidth = 0.000519861 - 8.58874e-06 * pt

    # select D mesons mass window
    if selection == 'sgn':
        lower_selection = mass_mean - nsigma_mass * massWidth
        upper_selection = mass_mean + nsigma_mass * massWidth
    elif selection == 'sbl':
        lower_selection = mass_mean - nsigma_offset * massWidth - sideband_width
        upper_selection = mass_mean - nsigma_offset * massWidth
    elif selection == 'sbr':
        lower_selection = mass_mean + nsigma_offset * massWidth
        upper_selection = mass_mean + nsigma_offset * massWidth + sideband_width

        if pdg == 411 and mass > lower_Dstar_removal and mass < upper_Dstar_removal:
            return False

    if mass > lower_selection and mass < upper_selection:
        return True

    return False
