import sys
import pandas as pd
import argparse

from ROOT import TCanvas, TFile, gInterpreter, TF1, TGraphAsymmErrors, TDatabasePDG, kAzure, kGray, TLegend, gROOT, TLatex, TGraphErrors, EColor, TLine
gInterpreter.ProcessLine('#include "../../fempy/utils/functions.h"')
from ROOT import GeneralCoulombLednickyDoubleSource


def EvalError(graph, x):
    for iPoint in range(graph.GetN() -1):
        xLeft = graph.GetPointX(iPoint)
        xRight = graph.GetPointX(iPoint+1)

        if xLeft < x < xRight:
            yUncLeft = graph.GetErrorY(iPoint)
            yUncRight = graph.GetErrorY(iPoint + 1)

            m = (yUncRight - yUncLeft) / (xRight - xLeft)
            return m * (x - xLeft) + yUncLeft

    # return 100% error in case the point is outside the range of the graph
    print("x outside of the graph range")
    return graph.Eval(x)



def ComputePulls(gStat, gSyst, gTheo, plotRangeX):
    gPulls = TGraphErrors(1)
    iPull = 0
    for iPoint in range(gStat.GetN()):
        x = gStat.GetPointX(iPoint)
        if plotRangeX[0] < x < plotRangeX[1]:
            y = gStat.GetPointY(iPoint)
            yStat = gStat.GetErrorY(iPoint)
            ySyst = gSyst.GetErrorY(iPoint)
            yTheo = gTheo.Eval(x)
            yTheoUnc = EvalError(gTheo, x)
            print(y - yTheo - y)
            print('pull iPoint: ', iPoint, (y - yTheo)/(yTheoUnc ** 2 + yStat ** 2 + ySyst ** 2) ** 0.5, yTheoUnc)
            gPulls.SetPoint(iPull, x, (y - yTheo)/(yTheoUnc ** 2 + yStat ** 2 + ySyst ** 2) ** 0.5)
            iPull += 1
    return gPulls


# draw yuki yuki
def ComputeChi2(gStat, gSyst, gTheo, plotRangeX):
    chi2 = 0
    ndf = 0
    for iPoint in range(gStat.GetN()):
        x = gStat.GetPointX(iPoint)
        if plotRangeX[0] < x and x < plotRangeX[1]:
            y = gStat.GetPointY(iPoint)
            yStat = gStat.GetErrorY(iPoint)
            ySyst = gSyst.GetErrorY(iPoint)
            yTheo = gTheo.Eval(x)
            yTheoUnc = EvalError(gTheo, x)
            print('chi2 iPoint: ', iPoint, (y - yTheo)/(yTheoUnc ** 2 + yStat ** 2 + ySyst ** 2) ** 0.5)
            chi2 += (y - yTheo) ** 2 / (yTheoUnc ** 2 + yStat ** 2 + ySyst ** 2)
            ndf += 1

    return chi2, ndf


def ComputeSystUnc(gTot, gStat):
    if gTot.GetN() != gStat.GetN():
        print("Error. different number of points. Exit!")
        sys.exit()
    gSyst = TGraphErrors(1)
    
    for iPoint in range(gTot.GetN()):
        xTot = gTot.GetPointX(iPoint)
        xStat = gStat.GetPointX(iPoint)
        if (abs((xStat - xTot)/xTot) > 1.e-6):
            print('x values are different for stat and syst unc graphs')

        
        yTot = gTot.GetPointY(iPoint)
        yStat = gStat.GetPointY(iPoint)
        # if (abs((yStat - yTot)/yTot) > 1.e-6):
        #     print('y values are different for stat and syst unc graphs')
        
        
        xUnc = gStat.GetErrorX(iPoint)
        yTotUnc = gTot.GetErrorY(iPoint)
        yStatUnc = gStat.GetErrorY(iPoint)
        shift = abs(yTot - yStat)
        if yTotUnc > yStatUnc:
            ySystUnc = (shift **2 + yTotUnc ** 2 - yStatUnc ** 2) ** 0.5
        else:
            print("Warning: stat error larger than total error! setting syst error = y_shift.")
            ySystUnc = shift
        gSyst.SetPoint(iPoint, xStat, yStat)
        gSyst.SetPointError(iPoint, xUnc, ySystUnc)

    return gSyst



gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument('pair')
parser.add_argument('comb')
args = parser.parse_args()

plotRangeX = [0, 500]
plotRangeY = [0.6, 1.4]

cCF = TCanvas('cCF', '', 600, 900)
cCF.Divide(1, 2)
pad = cCF.cd(1)
pad.DrawFrame(plotRangeX[0], plotRangeY[0], plotRangeX[1], plotRangeY[1], ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')


# draw Lednicky
fLL = TF1('fLL', GeneralCoulombLednickyDoubleSource, 0.1, 300, 8)
if args.pair == 'DPi':
    gYuki = TGraphAsymmErrors("yuki/corr_model_unc/corr_unc_Liu/PipDp.dat", "%lg %lg %lg %lg")
    w1 = 0.66
    s1 = 0.97
    s2 = 2.52

    gStat = TFile("/home/daniel/paper/CharmPaper/figures/final_D/Result_PipDp.root").Get('CF_corr_NBL_truth')
    gTot = TFile("/home/daniel/paper/CharmPaper/figures/final_D/Result_PipDp_tot.root").Get('CF_corr_NBL_truth')
    gSyst = ComputeSystUnc(gTot, gStat)


    
elif args.pair == 'DK':
    gYuki = TGraphAsymmErrors("yuki/corr_model_unc/corr_unc_Liu/KpDp.dat", "%lg %lg %lg %lg")
    w1 = 0.78
    s1 = 0.86
    s2 = 2.03
    
gYuki.SetFillColorAlpha(kAzure+2, 0.5)
gYuki.SetLineColor(1)
gYuki.SetMarkerColor(1)
gYuki.Draw('same ce3')
    
a0 = -0.06
m1 = TDatabasePDG.Instance().GetParticle(411).Mass()
m2 = TDatabasePDG.Instance().GetParticle(211 if args.pair == 'DPi' else 321).Mass()
m = m1 * m2 / (m1 + m2) * 1000

fLL.SetParameter(0, s1)
fLL.SetParameter(1, s2)
fLL.SetParameter(2, w1)
fLL.SetParameter(3, a0)
fLL.SetParameter(4, 0)
fLL.SetParameter(5, 0)
fLL.SetParameter(6, m)
fLL.SetParameter(7, +1)
fLL.Draw('same')

gSyst.SetFillColor(kGray)
gSyst.Draw('e2 same')
gStat.Draw('same')


leg = TLegend(0.5, 0.6, 0.8, 0.8)
if args.pair + args.comb == 'DPisc':
    header = 'D^{+}#pi^{+} #oplus D^{-}#pi^{-}'
elif  args.pair + args.comb == 'DKsc':
    header = 'D^{+}K^{+} #oplus D^{-}K^{-}'

leg.SetHeader(header)
leg.AddEntry(gStat, 'Data')
leg.AddEntry(gYuki, 'yuki (Liu) a0 = -0.1')

leg.AddEntry(fLL, f'Lednicky a_{{0}} = {a0}')
leg.Draw('same')

tl = TLatex()
chi2, ndf = ComputeChi2(gStat, gSyst, gYuki, plotRangeX)

tl.DrawLatexNDC(0.15, 0.8,f'#chi^{{2}}/ndf = {chi2:.2f}/{ndf}')

line = TLine(plotRangeX[0], 1, plotRangeX[1], 1)
line.SetLineColor(EColor.kGray+2)
line.SetLineStyle(9)
line.Draw()

pad = cCF.cd(2)
pad.DrawFrame(plotRangeX[0], -4, plotRangeX[1], 4, ';#it{k}* (MeV/#it{c});Pull')
pad.SetGridy(True)

gPulls = ComputePulls(gStat, gSyst, gYuki, plotRangeX)
gPulls.SetMarkerColor(EColor.kBlack)
gPulls.SetLineColor(EColor.kBlack)
gPulls.SetMarkerStyle(20)

gPulls.Draw('pe same')



cCF.SaveAs('canvas.png')
cCF.SaveAs('canvas.pdf')

# Draw ratio LL/yuki
gRatioLLOverYuki = TGraphErrors(1)
gRatioLLOverYukiInteraction = TGraphErrors(1)
gLLMinusYuki = TGraphErrors(1)
for iPoint in range(gYuki.GetN()):
    x = gYuki.GetPointX(iPoint)
    yYuki = gYuki.GetPointY(iPoint)

    yLL = fLL.Eval(x)
    gRatioLLOverYuki.SetPoint(iPoint, x, yLL/yYuki)
    gRatioLLOverYukiInteraction.SetPoint(iPoint, x, (yLL-1)/(yYuki-1))
    gLLMinusYuki.SetPoint(iPoint, x, yLL-yYuki)
cRatioLLOverYuki = TCanvas('cRatioLLOverYuki', '', 600, 600)
cRatioLLOverYuki.SetLeftMargin(0.15)
gRatioLLOverYuki.SetLineWidth(2)
gRatioLLOverYuki.SetLineColor(EColor.kRed+2)
gRatioLLOverYuki.Draw('al')
gRatioLLOverYuki.SetTitle('LL vs Yuki;#it{k}* (MeV/#it{c});#it{C}_{LL}/#it{C}_{Yuki}')
cRatioLLOverYuki.SaveAs('RatioLLOverYuki.pdf')

cRatioLLOverYukiInteraction = TCanvas('cRatioLLOverYukiInteraction', '', 600, 600)
cRatioLLOverYukiInteraction.SetLeftMargin(0.15)
gRatioLLOverYukiInteraction.SetLineWidth(2)
gRatioLLOverYukiInteraction.SetLineColor(EColor.kRed+2)
gRatioLLOverYukiInteraction.SetTitle('LL vs Yuki;#it{k}* (MeV/#it{c});(#it{C}_{LL}-1)/(#it{C}_{Yuki}-1)')
gRatioLLOverYukiInteraction.Draw('al')
cRatioLLOverYukiInteraction.SaveAs('cRatioLLOverYukiInteraction.pdf')

cLLMinusYuki = TCanvas('cLLMinusYuki', '', 600, 600)
cLLMinusYuki.SetLeftMargin(0.15)
gLLMinusYuki.SetLineWidth(2)
gLLMinusYuki.SetLineColor(EColor.kRed+2)
gLLMinusYuki.SetTitle('LL vs Yuki;#it{k}* (MeV/#it{c});#it{C}_{LL} - #it{C}_{Yuki}')
gLLMinusYuki.Draw('al')
cLLMinusYuki.SaveAs('cLLMinusYuki.pdf')
