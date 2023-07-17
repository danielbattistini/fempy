import argparse

import ctypes

import numpy as np

from ROOT import gInterpreter, TDatabasePDG, TFile, TCanvas, TLatex, gStyle, kBlue, kGray, TF1, kOrange, TLegend, TGraphErrors, kWhite, kRed, kGreen, kBlack, RDataFrame
gInterpreter.ProcessLine('#include "combfit/functions.h"')
from ROOT import GeneralCoulombLednickyTwoRadii, GeneralCoulombLednickySecondTwoRadii

from fempy import TranslateToLatex
from fempy.utils.io import GetGraphsInDir


def Setstyle():
    gStyle.SetTextFont(42)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)

    gStyle.Reset("Plain")
    gStyle.SetOptTitle(False)
    gStyle.SetTitleBorderSize(0)
    gStyle.SetOptStat(0)
    gStyle.SetCanvasColor(10)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetFrameLineWidth(1)
    gStyle.SetFrameFillColor(kWhite)
    gStyle.SetPadColor(10)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetHistLineWidth(1)
    gStyle.SetHistLineColor(kRed)
    gStyle.SetFuncWidth(2)
    gStyle.SetFuncColor(kGreen)
    gStyle.SetLineWidth(2)
    gStyle.SetLabelSize(0.045, "xyz")
    gStyle.SetLabelOffset(0.01, "y")
    gStyle.SetLabelOffset(0.01, "x")
    gStyle.SetLabelColor(kBlack, "xyz")
    gStyle.SetTitleSize(0.05, "xyz")
    gStyle.SetTitleOffset(1.25, "y")
    gStyle.SetTitleOffset(1.2, "x")
    gStyle.SetTitleFillColor(kWhite)
    gStyle.SetTextSizePixels(26)
    gStyle.SetTextFont(42)
    gStyle.SetLegendFillColor(kWhite)
    gStyle.SetLegendFont(42)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetErrorX(0.005)


def SetHistStyle(hist):
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelOffset(0.01)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelOffset(0.01)
    hist.GetYaxis().SetTitleOffset(1.25)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)

    
def GraphAverage(graphsStat, graphsTot, shift=False):
    nPoints = graphsTot[0].GetN()

    # Stat
    ysStat = [[] for _ in range(nPoints)]
    for graph in graphsStat:
        for iPoint in range(nPoints):
            ysStat[iPoint].append(graph.GetPointY(iPoint))
    yavgStat = [np.mean(ysStat[iPoint]) for iPoint in range(nPoints)]
    ystdStat = [np.std(ysStat[iPoint]) for iPoint in range(nPoints)]

    # Total
    ysTot = [[] for _ in range(nPoints)]
    for graph in graphsTot:
        for iPoint in range(nPoints):
            ysTot[iPoint].append(graph.GetPointY(iPoint))
    yavgTot = [np.mean(ysTot[iPoint]) for iPoint in range(nPoints)]
    ystdTot = [np.std(ysTot[iPoint]) for iPoint in range(nPoints)]

    # compute syst unc
    if shift:
        shifts = [abs(tot - stat) for (stat, tot) in zip(yavgStat, yavgTot)]
        ystdSyst = [np.sqrt(tot**2 - stat**2 + s**2) for (stat, tot, s) in zip(ystdStat, ystdTot, shifts)]
    else:
        ystdSyst = [np.sqrt(tot**2 - stat**2) for (stat, tot) in zip(ystdStat, ystdTot)]

    gAverageStat = TGraphErrors(1)
    gAverageSyst = TGraphErrors(1)
    for iPoint in range(nPoints):
        gAverageStat.SetPoint(iPoint, graphsTot[0].GetPointX(iPoint), yavgStat[iPoint])
        gAverageSyst.SetPoint(iPoint, graphsTot[0].GetPointX(iPoint), yavgStat[iPoint])

        gAverageStat.SetPointError(iPoint, graphsTot[0].GetErrorX(iPoint), ystdStat[iPoint])
        gAverageSyst.SetPointError(iPoint, 0.5 * graphsTot[0].GetErrorX(iPoint), ystdSyst[iPoint])
    return {
        'stat': gAverageStat,
        'syst': gAverageSyst,
    }

def FunctionAverage(functionsStat, functionsTot, nPoints = 500, shift=False):
    yStat = [[] for _ in range(nPoints)]
    ySyst = [[] for _ in range(nPoints)]
    yTot = [[] for _ in range(nPoints)]

    xmin = ctypes.c_double()
    xmax = ctypes.c_double()
    functionsStat[0].GetRange(xmin, xmax)
    xx = np.linspace(xmin.value, xmax.value, nPoints)
    
    for func in functionsStat:
        for iPoint, x in enumerate(xx):
            yStat[iPoint].append(func.Eval(x))

    for func in functionsTot:
        for iPoint, x in enumerate(xx):
            yTot[iPoint].append(func.Eval(x))
    
    yavgStat = [np.mean(yStat[iPoint]) for iPoint in range(nPoints)]
    ystdStat = [np.std(yStat[iPoint]) for iPoint in range(nPoints)]
    yavgTot = [np.mean(yTot[iPoint]) for iPoint in range(nPoints)]
    ystdTot = [np.std(yTot[iPoint]) for iPoint in range(nPoints)]

    shifts = [abs(tot - stat) for tot, stat in zip(yavgTot, yavgStat)]
    ystdSystPlusShift = [np.sqrt(tot**2 + s**2) for tot, stat, s in zip(ystdTot, ystdStat, shifts)]
    ystd = ystdSystPlusShift if shift else ystdTot
    gAverage = TGraphErrors(1)
    
    for iPoint, (x, y, yerr) in enumerate(zip(xx, yavgStat, ystd)):
        gAverage.SetPoint(iPoint, x, y)
        gAverage.SetPointError(iPoint, 0, yerr)

    return gAverage


def ComputeBinBrackets(hist):
    gBrackets = TGraphErrors(1)
    for iBin in range(hist.GetNbinsX()):
        gBrackets.SetPoint(iBin, hist.GetBinCenter(iBin+1), hist.GetBinContent(iBin+1))
        gBrackets.SetPointError(iBin, hist.GetBinWidth(iBin+1)/2, 0)
    return gBrackets


def WeightedCoulombLednicky(x, par):
    # effective source radii
    r1 = par[0]
    r2 = par[1]

    # source weights
    w1 = par[2]
    w2 = 1 - w1

    scattLen = par[3]
    effRange = par[4]

    quantumStat = bool(par[5])
    redMass = par[6]
    chargeProd = par[7]

    # coulomb + strong CF
    gcl1 = GeneralCoulombLednicky(x[0], r1, scattLen, effRange, quantumStat, redMass, chargeProd)
    gcl2 = GeneralCoulombLednicky(x[0], r2, scattLen, effRange, quantumStat, redMass, chargeProd)

    return w1 * gcl1 + w2 * gcl2


def LoadRealCoulomb():
    fileSC = TFile('~/alice/CharmingAnalyses/DKDpi/feeddownContribCF/momRes/Dstar_PiplusDplusOutput.root')
    gSC = fileSC.Get('genuineCF')
    gSC.SetName('gCoulombSC')
    
    fileOC = TFile('~/alice/CharmingAnalyses/DKDpi/feeddownContribCF/momRes/Dstar_PiplusDminusOutput.root')
    gOC = fileOC.Get('genuineCF')
    gOC.SetName('gCoulombOC')
    
    return {
        'sc': gSC,
        'oc': gOC,
    }

def RedMass(m1, m2):
    return m1 * m2 / (m1 + m2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile', nargs='?', default='/home/daniel/an/DstarPi/20_luuksel/SimFit_nopc_kStarBW50MeV_bs10002syst_chi2lt10_radii.root')
    args = parser.parse_args()

    lightMass = TDatabasePDG.Instance().GetParticle(211).Mass()
    heavyMass = TDatabasePDG.Instance().GetParticle(413).Mass()

    inFile = TFile(args.inFile)

    graphsTot = GetGraphsInDir(inFile.Get('tot/iters'))
    graphsSCTot = [g for g in graphsTot if 'SC' in g.GetName()]
    graphsOCTot = [g for g in graphsTot if 'OC' in g.GetName()]

    graphs = GetGraphsInDir(inFile.Get('stat/iters'))
    graphsSCStat = [g for g in graphs if 'SC' in g.GetName()]
    graphsOCStat = [g for g in graphs if 'OC' in g.GetName()]
    
    dgCF = {
        'sc': GraphAverage(graphsSCStat, graphsSCTot, shift=True),
        'oc': GraphAverage(graphsOCStat, graphsOCTot, shift=True),
    }

    # build LL curves
    paramsLLTot = np.array(list(RDataFrame(inFile.Get('tot/tResults')).AsNumpy(['a0sin', 'a0tri', 'r1', 'r2', 'w1']).values())).T
    paramsLLStat = np.array(list(RDataFrame(inFile.Get('stat/tResults')).AsNumpy(['a0sin', 'a0tri', 'r1', 'r2', 'w1']).values())).T

    lfLLSCStat = []
    for a0sinStat, a0triStat, r1, r2, w1 in paramsLLStat:
        lfLLSCStat.append(TF1("fLLsc", GeneralCoulombLednickyTwoRadii, 4, 300, 8))
        lfLLSCStat[-1].FixParameter(0,r1)
        lfLLSCStat[-1].FixParameter(1,r2)
        lfLLSCStat[-1].FixParameter(2,w1)
        lfLLSCStat[-1].FixParameter(3,a0sinStat)
        lfLLSCStat[-1].FixParameter(4,0)
        lfLLSCStat[-1].FixParameter(5,0)
        lfLLSCStat[-1].FixParameter(6,RedMass(lightMass, heavyMass)*1000)
        lfLLSCStat[-1].FixParameter(7,1)

    lfLLOCStat = []
    for a0sinStat, a0triStat, r1, r2, w1 in paramsLLStat:
        lfLLOCStat.append(TF1("fLLsc", GeneralCoulombLednickySecondTwoRadii, 4, 300, 10))
        lfLLOCStat[-1].FixParameter(0,r1)
        lfLLOCStat[-1].FixParameter(1,r2)
        lfLLOCStat[-1].FixParameter(2,w1)
        lfLLOCStat[-1].FixParameter(3,a0sinStat)
        lfLLOCStat[-1].FixParameter(4,0)
        lfLLOCStat[-1].FixParameter(5,a0triStat)
        lfLLOCStat[-1].FixParameter(6,0)
        lfLLOCStat[-1].FixParameter(7,0)
        lfLLOCStat[-1].FixParameter(8,RedMass(lightMass, heavyMass)*1000)
        lfLLOCStat[-1].FixParameter(9,-1)
    
    lfLLSCTot = []
    for a0sinTot, a0triTot, r1, r2, w1 in paramsLLTot:
        lfLLSCTot.append(TF1("fLLsc", GeneralCoulombLednickyTwoRadii, 4, 300, 8))
        lfLLSCTot[-1].FixParameter(0,r1)
        lfLLSCTot[-1].FixParameter(1,r2)
        lfLLSCTot[-1].FixParameter(2,w1)
        lfLLSCTot[-1].FixParameter(3,a0sinTot)
        lfLLSCTot[-1].FixParameter(4,0)
        lfLLSCTot[-1].FixParameter(5,0)
        lfLLSCTot[-1].FixParameter(6,RedMass(lightMass, heavyMass)*1000)
        lfLLSCTot[-1].FixParameter(7,1)

    lfLLOCTot = []
    for a0sinTot, a0triTot, r1, r2, w1 in paramsLLTot:
        lfLLOCTot.append(TF1("fLLsc", GeneralCoulombLednickySecondTwoRadii, 4, 300, 10))
        lfLLOCTot[-1].FixParameter(0,r1)
        lfLLOCTot[-1].FixParameter(1,r2)
        lfLLOCTot[-1].FixParameter(2,w1)
        lfLLOCTot[-1].FixParameter(3,a0sinTot)
        lfLLOCTot[-1].FixParameter(4,0)
        lfLLOCTot[-1].FixParameter(5,a0triTot)
        lfLLOCTot[-1].FixParameter(6,0)
        lfLLOCTot[-1].FixParameter(7,0)
        lfLLOCTot[-1].FixParameter(8,RedMass(lightMass, heavyMass)*1000)
        lfLLOCTot[-1].FixParameter(9,-1)
        
        
                    
    dgLL = {
        'sc': FunctionAverage(lfLLSCStat, lfLLSCTot, shift=True),
        'oc': FunctionAverage(lfLLOCStat, lfLLOCTot, shift=True),
    }
    _, _, r1, r2, w1 = paramsLLTot[0]

    a0sinTot = np.mean(paramsLLTot.T[0])
    a0sinStat = np.mean(paramsLLStat.T[0])
    a0sinUncStat = np.std(paramsLLStat.T[0])
    a0sinUncSyst = np.sqrt(np.std(paramsLLTot.T[0])**2 - a0sinUncStat**2)

    a0triTot = np.mean(paramsLLTot.T[1])
    a0triStat = np.mean(paramsLLStat.T[1])
    a0triUncStat = np.std(paramsLLStat.T[1])
    a0triUncSyst = np.sqrt(np.std(paramsLLTot.T[1])**2 - a0triUncStat**2)

    print(f'a0(3/2) = {a0sinStat:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    print(f'a0(1/2) = {a0triStat:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')

    dgLLCoulombOnly = LoadRealCoulomb()

    oFileNameBase = '/home/daniel/an/DstarPi/20_luuksel/SimFitPlot_DstarPi_simfit_shift'
    oFile = TFile(f'{oFileNameBase}.root', 'create')

    Setstyle()
    cCF = TCanvas('cCF', '', 1000, 500)
    cCF.Divide(2, 1)
    for iPad, comb in enumerate(['sc', 'oc']):
        oFile.mkdir(comb)
        oFile.cd(comb)

        pad = cCF.cd(iPad+1)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.03)
        pad.DrawFrame(0, 0.7, 299.9999, 1.8999, TranslateToLatex(';__kStarMeV__;__C__'))

        # Draw LL
        dgLL[comb].SetFillColorAlpha(kBlue+1, 0.7)
        dgLL[comb].SetLineColorAlpha(0, 0)
        dgLL[comb].Draw('same e3')
        
        dgLLCoulombOnly[comb].SetLineColor(kOrange+7)
        dgLLCoulombOnly[comb].SetLineStyle(1)
        dgLLCoulombOnly[comb].SetLineWidth(2)
        dgLLCoulombOnly[comb].Draw('same l')

        # Draw syst
        dgCF[comb]['syst'].SetFillColorAlpha(kGray+2, 0.65)
        dgCF[comb]['syst'].Draw('same e2')
        dgCF[comb]['syst'].Write('gCFSyst')
        # Draw stat
        dgCF[comb]['stat'].SetMarkerSize(1)
        dgCF[comb]['stat'].SetMarkerStyle(24)
        dgCF[comb]['stat'].Draw('same pez')
        dgCF[comb]['stat'].Write('gCFStat')

        # brackets
        gBrackets = TGraphErrors(1)
        for iPoint in range(6):
            gBrackets.SetPoint(iPoint, 50 * iPoint + 25, dgCF[comb]['stat'].GetPointY(iPoint))
            gBrackets.SetPointError(iPoint, 25, 0)

        SetHistStyle(gBrackets)
        gBrackets.SetLineWidth(1)
        gBrackets.SetLineColorAlpha(kBlack, 0.9)
        gBrackets.DrawClone("same []")
        gBrackets.Write('gBrackets')

        tl = TLatex()
        tl.SetTextSize(0.04)
        tl.SetNDC(True)
        if iPad == 0:
            tl.DrawLatex(0.22, 0.88, 'ALICE pp #sqrt{#it{s}} = 13 TeV')
            tl.DrawLatex(0.22, 0.82, 'High-mult. (0 #minus 0.17% INEL > 0)')
            tl.DrawLatex(0.22, 0.76, TranslateToLatex(f'kDstarPi_{comb}'))
            
            leg = TLegend(0.2, 0.55, 0.9, 0.7)
            leg.SetTextSizePixels(12)
            leg.SetTextSize(0.035)
            leg.SetLineColorAlpha(0, 0)
            leg.SetFillStyle(0)
            leg.AddEntry(dgCF['sc']['stat'], 'Genuine CF', 'pel')
            leg.AddEntry(dgLLCoulombOnly['sc'], 'Coulomb only', 'l')
            leg.AddEntry(dgLL['sc'], 'Lednick#acute{y}-Lyuboshits model', 'f')
            leg.Draw()
        else:
            tl.DrawLatex(0.20, 0.88, TranslateToLatex(f'kDstarPi_{comb}'))
            tl.DrawLatex(0.20, 0.25, f'a_{{0}}(I=3/2) = {a0sinStat:.2f} #pm {a0sinUncStat:.2f} (stat) #pm {a0sinUncSyst:.2f} (syst) fm')
            tl.DrawLatex(0.20, 0.25-0.05, f'a_{{0}}(I=1/2) = {a0triStat:.2f} #pm {a0triUncStat:.2f} (stat) #pm {a0triUncSyst:.2f} (syst) fm')
        
    for ext in ['pdf', 'png', 'eps']:
        cCF.SaveAs(f'{oFileNameBase}.{ext}')
    oFile.cd('/')
    cCF.Write()