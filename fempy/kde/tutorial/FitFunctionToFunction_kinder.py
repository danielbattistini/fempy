import numpy as np
from ROOT import TF1, TF2, gRandom, TH1D, TCanvas, TKDE, kBlue, kGreen, gROOT, TLegend, gStyle, gInterpreter, gSystem, kViolet
from sys import exit

gInterpreter.ProcessLine('#include "../utils/FitFunctionToFunction.h"')
gSystem.Load('../utils/FitFunctionToFunction_cxx.so')

from ROOT import FitFunctionToFunction, variableSettings

fTrue = TF1('fTrue', 'x + x*x+10', 0, 10)

constunc = False
if constunc:
    fUncUpper = TF1('fUncUpper', 'fTrue + 5', 0, 10)
    fUncLower = TF1('fUncLower', 'fTrue - 5', 0, 10)
else:
    fUncUpper = TF1('fUncUpper', 'fTrue * 1.1', 0, 10)
    fUncLower = TF1('fUncLower', 'fTrue * 0.9', 0, 10)

fFit = TF1('fFit', '[0]+[1]*x', 0, 10)

nPar = 2
fitter = FitFunctionToFunction(fTrue, fFit, fUncUpper, fUncLower, nPar)

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

fitter.SetVariableSettings(0, p0)
fitter.SetVariableSettings(1, p1)
fitter.SetFitRange(0, 10)

fitter.Fit()
# [print(f.GetName()) for f in list(gROOT.GetListOfFunctions())]

# cd = TCanvas('cd', 'cd', 600, 600)
# fDelta = TF1('f', '4*((fFunc1-fFunc2)/(fFunc2UpperUnc - fFunc2LowerUnc))**2', 0, 10)
# # fDelta = TF1('f', '((fFunc1-fFunc2)/(fFunc2UpperUnc - fFunc2LowerUnc)**2', 0, 10)
# # fDelta = TF1('f', '((fTrue-fFit)/(fUncUpper - fUncLower))*((fTrue-fFit)/(fUncUpper - fUncLower))', 0, 10)
# fDelta.SetLineColor(kViolet)
# fDelta.Draw()

c = TCanvas('c', 'c', 600, 600)
fTrue.GetYaxis().SetRangeUser(-20, 130)
fTrue.Draw()

fUncUpper.SetLineStyle(2)
fUncUpper.Draw('same')
fUncLower.SetLineStyle(2)
fUncLower.Draw('same')

fFit.SetLineColor(kBlue)
fFit.Draw('same')

print('---------->>>> Chi2: ', fitter.GetChi2())

leg = TLegend(0.2, 0.65, 0.5, 0.85)
leg.AddEntry(fTrue, 'function to fit')
leg.AddEntry(fUncLower, 'uncertainty (1#sigma)')
leg.AddEntry(fFit, 'fit function (pol1)')
leg.Draw('same')

if constunc:
    c.SaveAs('fit_function2function_constunc.pdf')
else:
    c.SaveAs('fit_function2function_varunc.pdf')