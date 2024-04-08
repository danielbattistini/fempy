import pandas as pd
import numpy as np
import argparse

from ROOT import TFile, TCanvas, TLegend, EColor, TGraph, TH2D, gStyle
gStyle.SetPalette(55)

parser = argparse.ArgumentParser()
parser.add_argument('inFile')
parser.add_argument('oFileWoExt')
args = parser.parse_args()

inFile = TFile(args.inFile)
tResults = inFile.Get('tResults')

dfScatLenVsPotAtr = pd.read_csv('../theory/pot2a0/f0_Dpi_atr.dat', header=None, delim_whitespace=True)
# reverse the order of the dataframe in order to avoid root getting confused with  direction of connecting lines
dfScatLenVsPotAtr = dfScatLenVsPotAtr.iloc[::-1]

dfScatLenVsPotRep = pd.read_csv('../theory/pot2a0/f0_Dpi_rep.dat', header=None, delim_whitespace=True)
dfScatLenVsPot = pd.concat([dfScatLenVsPotAtr, dfScatLenVsPotRep])

gPot2ScatLen = TGraph(len(dfScatLenVsPot), dfScatLenVsPot.to_numpy()[:, 0], dfScatLenVsPot.to_numpy()[:, 1])

hScatLen = TH2D('hScatLen', ';a_{0}^{D#pi}(I=3/2) (fm);a_{0}^{D#pi}(I=1/2) (fm)', 50, -0.15, 0.15, 50, -0.15, 0.15)
for iter in tResults:
    print(iter.vq, iter.vd)

    hScatLen.Fill(gPot2ScatLen.Eval(iter.vq), gPot2ScatLen.Eval(iter.vd))



# Scattering lenght
cScatLenCorr = TCanvas('cScatLenCorr', '', 600, 600)
cScatLenCorr.SetLeftMargin(0.15)
hScatLen.Draw('colz')
cScatLenCorr.SaveAs(f'{args.oFileWoExt}_ScatLenCorr.pdf')


# Potential
cPotCorr = TCanvas('cPotCorr', '', 600, 600)
cPotCorr.SetLeftMargin(0.15)
frame = cPotCorr.DrawFrame(-1500, -1500, 1500, 1500, 'Potential;V_{q} (MeV);V_{d} (MeV);Counts')
tResults.Draw('vd:vq', '', 'same colz')
cPotCorr.SaveAs(f'{args.oFileWoExt}_PotCorr.pdf')


# Potential vs scatering length
cPotVsScatLen = TCanvas('cPotVsScatLen', '', 600, 600)
cPotVsScatLen.SetLeftMargin(0.15)
frame = cPotVsScatLen.DrawFrame(-5001, -1.5, 5001, 1.5, 'Potential vs scat len;V_{0} (MeV);a_{0}^{D#pi} (fm)')
gPot2ScatLen.SetLineWidth(2)
gPot2ScatLen.SetLineColor(EColor.kRed+2)

gPot2ScatLen.Draw('same l')
