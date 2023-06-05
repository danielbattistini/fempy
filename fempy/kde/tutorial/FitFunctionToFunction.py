from os import ftruncate
import numpy as np
from ROOT import TF1, TF2, gRandom, TH1D, TCanvas, TKDE, kBlue, kGreen, gROOT, TLegend, gStyle, gInterpreter, gSystem, SetOwnership
from sys import exit

gInterpreter.ProcessLine('#include "../utils/FitFunctionToFunction.h"')
gSystem.Load('../utils/FitFunctionToFunction_cxx.so')

from ROOT import FitFunctionToFunction, variableSettings
nMeas = 100

TF1.DefaultAddToGlobalList(True)

fTrue = TF1("fTrue", '[0]*pow(e, -x/[1])', 0, 10)

fTrue.SetParameter(0, 1)
fTrue.SetParameter(1, 2.5)

data = []
hData = TH1D('hData', 'hData', 50, 0, 10)
for _ in range(nMeas):
    xTmp = fTrue.GetRandom()
    data.append(xTmp)    
    hData.Fill(xTmp)

# make kde estimation
dataArray = np.array(data, 'd')
kde = TKDE(nMeas, dataArray, 0, 10, "KernelType:Gaussian;Iteration:Adaptive;Mirror:MirrorAsymBoth;Binning:RelaxedBinning")

fKde = kde.GetFunction()
fKde.AddToGlobalList()


fFit = TF1("fFit", '[0]*pow(e, -x/[1])', 0, 10)


fLine = TF1("fLine", "[0]+[1]*x", 0, 10)
# fKde = fTrue
# fFit = fLine
# fFit = fLine
# fFit.SetParLimits(0, 0, 1.5)
# fFit.SetParLimits(1, -10, -0.2)
# fFit.SetParameters(0, 1)
# fFit.SetParameters(1, 0.2)
# print(fKde.GetExpFormula())
nPar = 2
fUpperUnc = kde.GetUpperFunction(0.683)
fLowerUnc = kde.GetLowerFunction(0.683)

fUpperUnc.AddToGlobalList()
fLowerUnc.AddToGlobalList()

SetOwnership(fKde, False)
SetOwnership(fUpperUnc, False)
SetOwnership(fLowerUnc, False)

fFit = TF1('fFit', '[0]+[1]*x', 0, 10)
fTrue = TF1('fTrue', '1 -x + x*x', 0, 10)
constunc = True
# if constunc:
#     fUncUpper = TF1('fUncUpper', 'fTrue + 5', 0, 10)
#     fUncLower = TF1('fUncLower', 'fTrue - 5', 0, 10)
# else:
#     fUncUpper = TF1('fUncUpper', 'fTrue * 1.1', 0, 10)
#     fUncLower = TF1('fUncLower', 'fTrue * 0.9', 0, 10)
# fitter = FitFunctionToFunction(fTrue, fFit, fUncUpper, fUncLower, nPar)
fitter = FitFunctionToFunction(fTrue, fFit, fUpperUnc, fLowerUnc, nPar)


# fitter = FitFunctionToFunction(fKde, fFit, fUpperUnc, fLowerUnc, nPar)

p0 = variableSettings()
p0.fName = 'p0'
p0.fStart = 0
p0.fLimLow = -30
p0.fLimUp = 30
p0.fStep = 0.1

p1 = variableSettings()
p1.fName = 'p1'
p1.fStart = 1
p1.fLimLow = -30
p1.fLimUp = 20
p1.fStep = 0.1

# fitter = FitFunctionToFunction(fKde, fFit)

# fitter = FitFunctionToFunction(fKde, fFit)
# fitter = FitFunctionToFunction(fTrue, fLine)

fitter.SetVariableSettings(0, p0)
fitter.SetVariableSettings(1, p1)
fitter.SetFitRange(0, 10)

[print(f.GetName(), end='\t') for f in list(gROOT.GetListOfFunctions())]
print('\n\n')
[print(f.GetName(), end='\t') for f in list(gROOT.GetListOfGlobalFunctions())]

exit()
fitter.Fit()
# exit()
# print(f'p0: {fFit.GetParameter(0)}   p1: {fFit.GetParameter(1)}')
# exit()

# chi2 = fitter.GetChi2()
# print(chi2)

gStyle.SetOptStat(0)
c = TCanvas('c', 'c', 600, 600)

hData.Scale(fTrue.Integral(0, 10)/hData.Integral('width'))
hData.Draw('')
fTrue.Draw('same')

def scaledFunctionCentr(x, par):
    return fTrue.Integral(0, 10) * kde.GetFunction().Eval(x[0])

def scaledFunctionUpper(x, par):
    return fTrue.Integral(0, 10) * kde.GetUpperFunction(0.683).Eval(x[0])

def scaledFunctionLower(x, par):
    return fTrue.Integral(0, 10) * kde.GetLowerFunction(0.683).Eval(x[0])

scaledFuncCentr = TF1('cskjdc', scaledFunctionCentr, 0, 10)
scaledFuncUpper = TF1('cskjdc', scaledFunctionUpper, 0, 10)
scaledFuncLower = TF1('cskjdc', scaledFunctionLower, 0, 10)
scaledFuncCentr.SetLineColor(kBlue)
scaledFuncUpper.SetLineStyle(2)
scaledFuncLower.SetLineStyle(2)
scaledFuncUpper.SetLineColor(kBlue)
scaledFuncLower.SetLineColor(kBlue)
scaledFuncCentr.Draw('same')
scaledFuncUpper.Draw('same')
scaledFuncLower.Draw('same')

leg = TLegend(0.5, 0.7, 0.85, 0.85)
leg.AddEntry(hData, 'data')
leg.AddEntry(fTrue, 'true pdf')
leg.AddEntry(scaledFuncCentr, 'kde')
leg.AddEntry(fFit, 'fit to kde')
leg.Draw('same')


fFit.SetLineColor(kGreen)
fFit.Draw('same')
# [print(f.GetName()) for f in list(gROOT.GetListOfFunctions())]


# cc = TCanvas("c2", "c", 600, 600)

# print(fitter.customChi2Formula(1, 2))
# def customChi2form(xx, par):
#     return fitter.customChi2Formula(xx, par)

# fChi2_2D = TF2("fChi2", customChi2form, -2, 2, -2, 2)
# fChi2_2D.Draw("colz")

