'''
Format stings.
'''

import sys

from ROOT import gStyle, gROOT


def GetNormRangeFromLabel(pair):
    if 'DPi' in pair:
        return [1000, 1500]
    elif 'DK' in pair:
        return [1500, 2000]
    if 'DstarPi' in pair:
        return [1500, 2000]
    else:
        print("Error: pair not implemented. Exit!")
        sys.exit()


def GetNormRangeFromPDG(pdgHeavy, pdgLight):
    if abs(pdgHeavy) == 411 and abs(pdgLight) == 211:  # D pion
        return GetNormRangeFromLabel('DPi')
    elif abs(pdgHeavy) == 411 and abs(pdgLight) == 321:  # D kaon
        return GetNormRangeFromLabel('DK')
    elif abs(pdgHeavy) == 413 and abs(pdgLight) == 211:  # Dstar ion
        return GetNormRangeFromLabel('DstarPi')
    else:
        print("Error: pair not implemented. Exit!")
        sys.exit()


dPDG2Label = {
    411: 'D',
    413: 'Dstar',
    211: 'Pi',
    321: 'K',
}


def TranslateToLatex(text):
    '''
    The first charge always refers to the heavier particle.

    Charge combinations:
    pp: plus-plus
    mm: minus-minus
    pm: plus-minus
    mp: minus-plus
    sc: same-charge (pp+mm)
    oc: opposite-charge (pm+mp)
    '''

    dPairs2Latex = {
        # D-pion
        'kDPi_pp': "D^{+} #pi^{+}",
        'kDPi_mm': "D^{#minus} #pi^{#minus}",
        'kDPi_pm': "D^{+} #pi^{#minus}",
        'kDPi_mp': "D^{#minus} #pi^{+}",
        'kDPi_sc': "D^{+} #pi^{+} #oplus D^{#minus} #pi^{#minus}",
        'kDPi_oc': "D^{+} #pi^{#minus} #oplus D^{#minus} #pi^{+}",
        # D-kaon
        'kDK_pp': "D^{+} K^{+}",
        'kDK_mm': "D^{#minus} K^{#minus}",
        'kDK_pm': "D^{+} K^{#minus}",
        'kDK_mp': "D^{#minus} K^{+}",
        'kDK_sc': "D^{+} K^{+} #oplus D^{#minus} K^{#minus}",
        'kDK_oc': "D^{+} K^{#minus} #oplus D^{#minus} K^{+}",
        # Dstar-pion
        'kDstarPi_pp': "D*^{+} #pi^{+}",
        'kDstarPi_mm': "D*^{#minus} #pi^{#minus}",
        'kDstarPi_pm': "D*^{+} #pi^{#minus}",
        'kDstarPi_mp': "D*^{#minus} #pi^{+}",
        'kDstarPi_sc': "D*^{+} #pi^{+} #oplus D*^{#minus} #pi^{#minus}",
        'kDstarPi_oc': "D*^{+} #pi^{#minus} #oplus D*^{#minus} #pi^{+}",
        # axix labels
        'ciccio': "pasticcio",
        '__invmass_Dstar__': "#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})",
        '__invmass_D__': '#it{M}(K#pi#pi) (GeV/#it{c})^{2}',
        '__MeVc__': 'MeV/#it{c}',
        '__GeVc__': 'GeV/#it{c}',
        '__GeVc2__': 'GeV/#it{c}^{2}',
    }
    for key, value in dPairs2Latex.items():
        if key in text:
            text = text.replace(key, value)
    return text


def FigInit():
    # gROOT.SetBatch(True)

    gStyle.SetLineWidth(2)

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    gStyle.SetPadRightMargin(0.03)
    gStyle.SetPadLeftMargin(0.14)
    gStyle.SetPadTopMargin(0.03)
    gStyle.SetPadBottomMargin(0.1)
    
    gStyle.SetPalette(55)
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetLegendBorderSize(0)

    gStyle.SetTitleOffset(1., 'x')
    gStyle.SetTitleOffset(1.4, 'yz')
    
    gStyle.SetTitleSize(0.045, 'xyz')
    gStyle.SetLabelSize(0.04, 'xyz')

    gROOT.ForceStyle()
