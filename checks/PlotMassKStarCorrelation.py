from ROOT import TFile, TCanvas, gStyle

from fempy import TranslateToLatex

gStyle.SetPalette(55)

inFile = TFile('~/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
hhMassVsKStar = inFile.Get('pp/SE/hCharmMassVsKStar0')
cCorrelation = TCanvas('cCorr', '', 600, 600)
cCorrelation.SetLeftMargin(0.15)
cCorrelation.SetRightMargin(0.15)
cCorrelation.DrawFrame(0, 0.14, 200, 0.16, TranslateToLatex(';__kStarMeV__;__invmass_Dstar__;counts'))
hhMassVsKStar.Draw('same colz')
cCorrelation.SaveAs('~/an/DstarPi/20_luuksel/CharmMassVsKStar0_pp_SE.png')

