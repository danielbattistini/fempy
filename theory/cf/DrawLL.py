import argparse

from ROOT import TCanvas, gInterpreter, TF1, TDatabasePDG, TLatex, EColor, TLine

# gInterpreter.ProcessLine('#include "../../fempy/utils/functions.h"')
# from ROOT import GeneralCoulombLednickyDoubleSource

gInterpreter.ProcessLine('#include "../../combfit/functions.h"')
from ROOT import GeneralCoulombLednickyTwoRadii
from ROOT import ScalableGeneralCoulombLednickySecondTwoRadii

parser = argparse.ArgumentParser()
parser.add_argument('pair')
parser.add_argument('comb')
args = parser.parse_args()

plotRangeX = [0, 500]
plotRangeY = [0.6, 1.4]

# draw Lednicky
if args.pair == 'DPi':
    w1 = 0.66
    s1 = 0.97
    s2 = 2.52

elif args.pair == 'DK':
    w1 = 0.78
    s1 = 0.86
    s2 = 2.03

a0 = 0
m1 = TDatabasePDG.Instance().GetParticle(411).Mass()
m2 = TDatabasePDG.Instance().GetParticle(211 if args.pair == 'DPi' else 321).Mass()
m = m1 * m2 / (m1 + m2) * 1000

cCF = TCanvas('cCF', '', 600, 600)
cCF.DrawFrame(plotRangeX[0], plotRangeY[0], plotRangeX[1], plotRangeY[1], ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')

fLL = TF1('fLL', GeneralCoulombLednickyTwoRadii, 0.1, 300, 8)
fLL.SetParameter(0, s1)
fLL.SetParameter(1, s2)
fLL.SetParameter(2, w1)
fLL.SetParameter(3, a0)
fLL.SetParameter(4, 0)
fLL.SetParameter(5, 0)
fLL.SetParameter(6, m)
fLL.SetParameter(7, 1 if args.comb == 'sc' else -1)
fLL.Draw('same')

# fLL = TF1('fLL', GeneralCoulombLednickyDoubleSource, 0.1, 300, 8)
# fLL.SetParameter(0, s1)
# fLL.SetParameter(1, s2)
# fLL.SetParameter(2, w1)
# fLL.SetParameter(3, a0)
# fLL.SetParameter(4, 0)
# fLL.SetParameter(5, 0)
# fLL.SetParameter(6, m)
# fLL.SetParameter(7, 1 if args.comb == 'sc' else -1)
# fLL.Draw('same')

line = TLine(plotRangeX[0], 1, plotRangeX[1], 1)
line.SetLineColor(EColor.kGray+2)
line.SetLineStyle(9)
line.Draw()

fLLScaled = TF1('fLL', ScalableGeneralCoulombLednickySecondTwoRadii, 0.1, 300, 11)
fLLScaled.SetParameter(0, s1) # r1 (double)
fLLScaled.SetParameter(1, s2) # r2 (double)
fLLScaled.SetParameter(2, w1) # w1 (double)
fLLScaled.SetParameter(3, a0) # ScattLenSin (double)
fLLScaled.SetParameter(4, 0) # EffRangeSin (double)
fLLScaled.SetParameter(5, 0) # ScattLenTri (double)
fLLScaled.SetParameter(6, 0) # EffRangeTri (double)
fLLScaled.SetParameter(7, False) # QS (bool)
fLLScaled.SetParameter(8, m) # RedMass (double)
fLLScaled.SetParameter(9, 1 if args.comb == 'sc' else -1) # Q1Q (double)
fLLScaled.SetParameter(10, 0.99) # overall normalization (double)
fLLScaled.SetLineColor(EColor.kBlue)
fLLScaled.Draw("same")