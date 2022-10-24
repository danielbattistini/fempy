'''
Format stings.
'''

from ROOT import gStyle, gROOT

def TranslateTolatex(text):
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
    }
    for key in dPairs2Latex:
        if key in text:
            text = text.replace(key, dPairs2Latex[key])
    return text

def FigInit():
    gROOT.SetBatch(True)
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLineWidth(2)

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    gStyle.SetPadRightMargin(0.035)
    gStyle.SetPadLeftMargin(0.14)
    gStyle.SetPadTopMargin(0.035)
    gStyle.SetPadBottomMargin(0.1)
    gStyle.SetPalette(55)
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetLegendBorderSize(0)