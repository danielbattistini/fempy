from ROOT import TF1, TDatabasePDG, gInterpreter, TCanvas, TF1, kBlue, kBlue

gInterpreter.ProcessLine('#include "../fempy/utils/functions.h"')

from ROOT import GeneralCoulombLednicky, GeneralCoulombLednickyShort, TLegend, gROOT

mass_pi = TDatabasePDG.Instance().GetParticle(211).Mass()
mass_d = TDatabasePDG.Instance().GetParticle(411).Mass()
redmass = mass_pi*mass_d / (mass_pi + mass_d) * 1000

source = 1 # fm
a0 = +0.1 # fm

cLednicky = TCanvas('cLednicky', '', 600, 600)
cLednicky.DrawFrame(0, 0.8, 300, 1.2)
fGeneralCoulombLednickySC = TF1('fGeneralCoulombLednickySC', GeneralCoulombLednicky, 0.1, 300, 6)
fGeneralCoulombLednickySC.SetParameter(0, source) # source 1 fm
fGeneralCoulombLednickySC.SetParameter(1, a0)
fGeneralCoulombLednickySC.SetParameter(2, 0)
fGeneralCoulombLednickySC.SetParameter(3, 0)
fGeneralCoulombLednickySC.SetParameter(4, redmass)
fGeneralCoulombLednickySC.SetParameter(5, 1)
fGeneralCoulombLednickySC.Draw('same')

fGeneralCoulombLednickyShortSC = TF1('fGeneralCoulombLednickyShortSC', GeneralCoulombLednickyShort, 0.1, 300, 6)
fGeneralCoulombLednickyShortSC.SetParameter(0, source) # source 1 fm
fGeneralCoulombLednickyShortSC.SetParameter(1, a0)
fGeneralCoulombLednickyShortSC.SetParameter(2, 0)
fGeneralCoulombLednickyShortSC.SetParameter(3, 0)
fGeneralCoulombLednickyShortSC.SetParameter(4, redmass)
fGeneralCoulombLednickyShortSC.SetParameter(5, 1)
fGeneralCoulombLednickyShortSC.SetLineColor(kBlue)
fGeneralCoulombLednickyShortSC.Draw('same')

fGeneralCoulombLednickyOC = TF1('fGeneralCoulombLednickyOC', GeneralCoulombLednicky, 0.1, 300, 6)
fGeneralCoulombLednickyOC.SetParameter(0, source) # source 1 fm
fGeneralCoulombLednickyOC.SetParameter(1, a0)
fGeneralCoulombLednickyOC.SetParameter(2, 0)
fGeneralCoulombLednickyOC.SetParameter(3, 0)
fGeneralCoulombLednickyOC.SetParameter(4, redmass)
fGeneralCoulombLednickyOC.SetParameter(5, -1)
fGeneralCoulombLednickyOC.SetLineStyle(9)
fGeneralCoulombLednickyOC.Draw('same')

fGeneralCoulombLednickyShortOC = TF1('fGeneralCoulombLednickyShortOC', GeneralCoulombLednickyShort, 0.1, 300, 6)
fGeneralCoulombLednickyShortOC.SetParameter(0, source) # source 1 fm
fGeneralCoulombLednickyShortOC.SetParameter(1, a0)
fGeneralCoulombLednickyShortOC.SetParameter(2, 0)
fGeneralCoulombLednickyShortOC.SetParameter(3, 0)
fGeneralCoulombLednickyShortOC.SetParameter(4, redmass)
fGeneralCoulombLednickyShortOC.SetParameter(5, -1)
fGeneralCoulombLednickyShortOC.SetLineColor(kBlue)
fGeneralCoulombLednickyShortOC.SetLineStyle(9)
fGeneralCoulombLednickyShortOC.Draw('same')

leg = TLegend(0.4, 0.7, 0.9, 0.9)
leg.SetHeader(f'a_{{0}} = {a0}')
leg.AddEntry(fGeneralCoulombLednickyShortSC, 'LL basic same-charge')
leg.AddEntry(fGeneralCoulombLednickySC, 'LL additional terms same-charge')
leg.AddEntry(fGeneralCoulombLednickyShortOC, 'LL basic opposite-charge')
leg.AddEntry(fGeneralCoulombLednickyOC, 'LL additional terms opposite-charge')
leg.Draw('same')

cLednicky.SaveAs(f"cCompareLednickyFormulas_a0{a0}fm.png")