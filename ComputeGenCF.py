import sys
import yaml
import argparse
import copy

import pandas as pd
import numpy as np

from ROOT import TFile, TF1, TGraphErrors, TDatabasePDG, gInterpreter, TCanvas, kBlue, TLatex, gStyle, TLegend, gRandom, TNtuple, TH2D, RDataFrame, TH1D
gInterpreter.ProcessLine('#include "combfit/functions.h"')
gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter, ScalableGeneralCoulombLednickyTwoRadii, ScalableShifted300GeneralCoulombLednickyTwoRadii

LednickyModel = ScalableGeneralCoulombLednickyTwoRadii


import fempy
from fempy.utils.analysis import ComputeBinBrackets

gRandom.SetSeed(1)

def Average(hist, xmin, xmax):
    firstBin = hist.GetXaxis().FindBin(xmin * 1.0001)
    lastBin = hist.GetXaxis().FindBin(xmax * 0.9999)

    avg = 0
    for iBin in range(firstBin, lastBin + 1):
        avg += hist.GetBinContent(iBin) * hist.GetBinCenter(iBin)
    return avg/hist.Integral(firstBin, lastBin)


def StdDev(hist, xmin, xmax):
    firstBin = hist.GetXaxis().FindBin(xmin * 1.0001)
    lastBin = hist.GetXaxis().FindBin(xmax * 0.9999)
    mu = Average(hist, xmin, xmax)

    stdDev = 0
    for iBin in range(firstBin, lastBin + 1):
        stdDev += hist.GetBinContent(iBin) * (hist.GetBinCenter(iBin) - mu)**2
    return (stdDev/hist.Integral(firstBin, lastBin))**0.5


def ComputeNormFactor(se, me, start, end):
    firstBin = se.FindBin(start*1.0001)
    lastBin = me.FindBin(end*0.9999)
    return me.Integral(firstBin, lastBin) / se.Integral(firstBin, lastBin)


def SumLamPar(lam_par, treamtments):
    lam_par_summed = {
        'gen': 0,
        'flat': 0,
        'sb': 0,
    }
    for hk, h_lam in lam_par.items():
        for lk, l_lam in h_lam.items():
            lam_par_summed[treamtments[hk][lk]] += l_lam
    return lam_par_summed


def IsLamParMatValid(lam_par):
    total = 0
    for _, h_lam in lam_par.items():
        for _, l_lam in h_lam.items():
            total += l_lam
    return abs(total - 1) < 1e-6


def LoadLambdaParam(cfgCentr, npFracAbsVar=0., heavyPurity=None, lightFracRelVar=0):
    cfgVar = copy.deepcopy(cfgCentr)
    cfgVar['heavy'][1]['nonprompt']['frac'] += npFracAbsVar
    cfgVar['heavy'][0]['prompt']['frac'] = 1 - cfgVar['heavy'][1]['nonprompt']['frac']

    cfgVar['light'][1]['sec']['frac'] *= (1 + lightFracRelVar)
    cfgVar['light'][0]['prim']['frac'] = 1 - cfgVar['light'][1]['sec']['frac']

    if heavyPurity != None:
        cfgVar['heavy'][0]['prompt']['purity'] = heavyPurity
        cfgVar['heavy'][1]['nonprompt']['purity'] = heavyPurity
        cfgVar['heavy'][2]['bkg']['purity'] = 1. - heavyPurity

    lamParMatr = {}
    for heavyContrib in cfgVar['heavy']:
        heavyKey = list(heavyContrib.keys())[0]
        heavyPurity = heavyContrib[heavyKey]['purity']
        heavyFrac = heavyContrib[heavyKey]['frac']
        lamParMatr[heavyKey] = {}
        for lightContrib in cfgVar['light']:
            lightKey = list(lightContrib.keys())[0]
            lightPurity = lightContrib[lightKey]['purity']
            lightFrac = lightContrib[lightKey]['frac']

            lamParMatr[heavyKey][lightKey] = lightFrac * lightPurity * heavyFrac * heavyPurity
    return lamParMatr


def LoadGravities(hME, kStarBW, name='gGravities'):
    nKStarBins = round(hME.GetNbinsX()/kStarBW)
    xCF = [Average(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
    xCFUnc = [StdDev(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
    yCF = [1 for _ in range(nKStarBins)]
    yCFUnc = [0 for _ in range(nKStarBins)]
    gGravities = TGraphErrors(1)
    gGravities.SetName(name)
    for iBin, (x, y, xUnc, yUnc) in enumerate(zip(xCF, yCF, xCFUnc, yCFUnc)):
        gGravities.SetPoint(iBin, x, y)
        gGravities.SetPointError(iBin, xUnc, yUnc)
    return gGravities


def ApplyCenterOfGravity(hist, graph):
    nBins = hist.GetNbinsX()
    if nBins != graph.GetN():
        fempy.error(f"hist '{hist.GetName()}' has {hist.GetNbinsX()} bins "
                    f"but graph '{graph.GetName()}' has {graph.GetN()} points.")

    gCentered = graph.Clone(f'{hist.GetName()}_grav')
    for iBin in range(hist.GetNbinsX()):
        gCentered.SetPoint(iBin, graph.GetPointX(iBin), hist.GetBinContent(iBin+1))
        gCentered.SetPointError(iBin, graph.GetErrorX(iBin), hist.GetBinError(iBin+1))

    return gCentered


def RedMass(m1, m2):
    return m1 * m2 / (m1 + m2)


def MakeFlatHist(hSkeleton, value, name='hFlat'):
    hFlat = hSkeleton.Clone(name)
    hFlat.Reset()
    nBins = hFlat.GetNbinsX()

    # Set underflow and overflow to zero
    hFlat.SetBinContent(0, 0)
    hFlat.SetBinContent(nBins+1, 0)

    for iBin in range(nBins):
        hFlat.SetBinContent(iBin+1, value)
        hFlat.SetBinError(iBin+1, 0)

    return hFlat


def Bootstrap(hist):
    hBootstrapped = hist.Clone(f'{hist.GetName()}_bootstrap')
    for iBin in range(hist.GetNbinsX()+2):
        mu = hist.GetBinContent(iBin)
        sigma = hist.GetBinError(iBin)

        hBootstrapped.SetBinContent(iBin, gRandom.Gaus(mu, sigma))
        hBootstrapped.SetBinError(iBin, sigma)
    return hBootstrapped


def ComputePurity(hist, kStarBW, event, variation=""):
    nBins = round(3000/kStarBW)
    hPurity = TH1D(f'{hist.GetName()}_{event}_purity', '', nBins, 0, 3000)
    for iBin in range(nBins):
        firstBin = hist.FindBin(iBin * kStarBW * 1.0001)
        lastBin = hist.FindBin((iBin + 1) * kStarBW * 0.9999)

        hCharmMass = hist.ProjectionY('', firstBin, lastBin)
        fitter = MassFitter(hCharmMass, 'gaus', 'powex', 0.141, 0.154)
        fitter.Fit()

        nSigma = 2
        sgn = fitter.GetSignal(nSigma, 'data_minus_bkg')
        sgnUnc = fitter.GetSignalUnc(nSigma, 'data_minus_bkg')

        bkg = fitter.GetBackground(nSigma)
        bkgUnc = fitter.GetBackgroundUnc(nSigma)

        purity = sgn / (sgn + bkg)
        purityUnc = np.sqrt(bkg**2 * sgnUnc**2 + sgn**2 * bkgUnc ** 2) / (sgn + bkg)**2
        if variation == "":
            hPurity.SetBinContent(iBin+1, purity)
        elif variation == "+1":
            hPurity.SetBinContent(iBin+1, purity + purityUnc)
        elif variation == "-1":
            hPurity.SetBinContent(iBin+1, purity - purityUnc)

        hPurity.SetBinError(iBin+1, purityUnc)
    return hPurity


def VaryHistogram(hist, method):
    if method == 'bs':
        return Bootstrap(hist)
    elif isinstance(method, int):
        hVaried = hist.Clone()
        hVaried.Reset()
        for iBin in range(hVaried.GetNbinsX()+1):
            hVaried.SetBinContent(iBin, hist.GetBinContent(iBin) + method * hist.GetBinError(iBin))
            hVaried.SetBinError(iBin, hist.GetBinError(iBin))
        return hVaried
    else:
        fempy.error("not implemented")


def ComputeIntegratedPurity(hMass, fitRange=[0.141, 0.154]):
    firstBin = hMass.GetXaxis().FindBin(fitRange[0]*1.0001)
    lastBin = hMass.GetXaxis().FindBin(fitRange[1]*0.9999)

    particleYield = hMass.Integral(firstBin, lastBin)

    fitter = MassFitter(hMass, 'gaus', 'powex', fitRange[0], fitRange[1])
    fitter.Fit()

    sgn = fitter.GetSignal(2, 'data_minus_bkg')
    sgnUnc = fitter.GetSignalUnc(2, 'data_minus_bkg')

    bkg = fitter.GetBackground(2)

    purity = sgn / (sgn + bkg)
    purityUnc = purity * np.sqrt((sgnUnc/sgn)**2 + (1./np.sqrt(particleYield))**2)

    return (purity, purityUnc)


def ComputeScattPar(**kwargs):
    hCFSgn = kwargs['sgn']
    hCFSbr = kwargs['sbr']
    hCFMJ = kwargs['mj']
    lamPar = kwargs['lamPar']
    gGravities = kwargs['gGravities']
    scalableLL = kwargs['scalableLL']
    redMass = kwargs['redMass']
    radius1 = kwargs['radius1']
    radius2 = kwargs['radius2']
    weight1 = kwargs['weights1']
    fitRange = kwargs['fitRange']
    comb = kwargs['comb']
    tTrials = kwargs['tTrials']
    hhCFGen = kwargs['hhCFGen']
    iVar = kwargs['iVar']
    iIter = kwargs['iIter']
    purityVar = kwargs['purityVar']
    bkgFitRange = kwargs['bkgFitRange']

    if hCFSbr == None:
        # Compute the normalization of the MJ
        hCFNorm = hCFSgn / hCFMJ
        hCFNorm.SetName(f'hCFNorm{iIter}')
        hCFNorm.Write()

        # Fit the baseline
        fBaseLine = TF1(f'fBaseLine{iIter}', '[0]', 0, 3000)
        ApplyCenterOfGravity(hCFNorm, gGravities).Fit(fBaseLine, 'Q', '', bkgFitRange[0], bkgFitRange[1])
        blNorm = fBaseLine.GetParameter(0)

        # Compute the total background mode
        hCFBkg = blNorm * hCFMJ
        hCFBkg.SetName(f'hCFBkg{iIter}')
        hCFBkg.Write()

        # Compute the Gen CF
        hCFFlat = MakeFlatHist(hCFSgn, lamPar['flat'], 'hCFFlat')
        hCFGen = hCFSgn/(blNorm * hCFMJ) - hCFFlat
        hCFGen.Scale(1./lamPar['gen'])
        hCFGen.SetName(f'hCFGen{iIter}')
        hCFGen.Write()

    else:
        # Compute the normalization of the MJ
        hCFNorm = (hCFSgn - lamPar['sb'] * hCFSbr) / (hCFMJ * (lamPar['gen'] + lamPar['flat']))
        hCFNorm.SetName(f'hCFNorm{iIter}')
        hCFNorm.Write()

        # Fit the baseline
        fBaseLine = TF1(f'fBaseLine{iIter}', '[0]', 0, 3000)
        ApplyCenterOfGravity(hCFNorm, gGravities).Fit(fBaseLine, 'Q', '', bkgFitRange[0], bkgFitRange[1])
        blNorm = fBaseLine.GetParameter(0)

        # Compute the total background mode
        hCFBkg = lamPar['sb'] * hCFSbr + blNorm * hCFMJ * (lamPar['gen'] + lamPar['flat'])
        hCFBkg.SetName(f'hCFBkg{iIter}')
        hCFBkg.Write()

        # Compute the Gen CF
        hCFFlat = MakeFlatHist(hCFSgn, lamPar['flat'], 'hCFFlat')
        hCFGen = (hCFSgn - lamPar['sb'] * hCFSbr)/(blNorm * hCFMJ) - hCFFlat
        hCFGen.Scale(1./lamPar['gen'])
        hCFGen.SetName(f'hCFGen{iIter}')
        hCFGen.Write()

    # Add the CF to the hist2D
    if iIter > 0:
        for iBin in range(hCFGen.GetNbinsX()):
            hhCFGen.Fill(hCFGen.GetBinCenter(iBin+1), hCFGen.GetBinContent(iBin+1))

    # Compute the coulomb-only CF with Lednicky
    fCoulomb = TF1("fCoulomb", LednickyModel, fitRange[0], fitRange[1], 9)
    fCoulomb.FixParameter(0, radius1)
    fCoulomb.FixParameter(1, radius2)
    fCoulomb.FixParameter(2, weight1)
    fCoulomb.FixParameter(3, 0)
    fCoulomb.FixParameter(4, 0.)
    fCoulomb.FixParameter(5, 0.)
    fCoulomb.FixParameter(6, redMass)
    fCoulomb.FixParameter(7, 1 if comb == 'sc' else -1)
    fCoulomb.FixParameter(8, 1)
    fCoulomb.SetLineColor(kBlue)

    # Fit the CF
    gCFGen = ApplyCenterOfGravity(hCFGen, gGravities)
    fWeightedLL = TF1(f"fWeightedLL{iIter}", LednickyModel, fitRange[0], fitRange[1], 9)
    fWeightedLL.FixParameter(0, radius1)
    fWeightedLL.FixParameter(1, radius2)
    fWeightedLL.FixParameter(2, weight1)
    fWeightedLL.SetParameter(3, 0.)
    fWeightedLL.SetParLimits(3, -1, 1)
    fWeightedLL.FixParameter(4, 0.)
    fWeightedLL.FixParameter(5, 0.)
    fWeightedLL.FixParameter(6, redMass)
    fWeightedLL.FixParameter(7, 1 if comb == 'sc' else -1)
    if scalableLL:
        fWeightedLL.SetParameter(8, 1)
        fWeightedLL.SetParLimits(8, 0.1, 1.1)
    else:
        fWeightedLL.FixParameter(8, 1)

    status = gCFGen.Fit(fWeightedLL, "SMRQ+0").Status()
    chi2ndf = fWeightedLL.GetChisquare()/fWeightedLL.GetNDF()
    scattLen = fWeightedLL.GetParameter(3)
    scattLenUnc = fWeightedLL.GetParError(3)
    normLL = fWeightedLL.GetParameter(8)

    # Draw canvas
    cFit = TCanvas(f'cFit_{iIter}', '', 600, 600)
    cFit.DrawFrame(0, 0, 500, 2)
    hCFGen.Draw('pe')
    fWeightedLL.Draw('same')
    fCoulomb.Draw('same')

    step = 0.05
    tl = TLatex()
    tl.SetNDC()
    tl.SetTextSize(0.035)
    tl.SetTextFont(42)
    tl.DrawLatex(0.2, 0.85 - 0 * step, fempy.utils.format.TranslateToLatex(f'k{args.pair}_{comb}   iVar: {iVar}'))
    tl.DrawLatex(0.2, 0.85 - 1 * step, f'a_{{0}} = {scattLen:.3f} #pm {scattLenUnc:.3f}')
    tl.DrawLatex(0.2, 0.85 - 2 * step, f'#chi^{{2}}/ndf = {chi2ndf:.2f}')

    leg = TLegend(0.5, 0.7, 0.85, 0.85)
    leg.AddEntry(hCFGen, 'Data')
    leg.AddEntry(fWeightedLL, 'Fit LL')
    leg.AddEntry(fCoulomb, 'Coulomb LL')
    leg.Draw()
    if chi2ndf > 3 or abs(scattLen) > 0.99:
        tl.DrawLatex(0.2, 0.85 - 3 * step, 'BAD FIT')
    else:
        tTrials.Fill(scattLen, scattLenUnc, status, chi2ndf, iVar, purityVar, weight1, radius1, radius2, fitRange[1], bkgFitRange[0], bkgFitRange[1], lamPar['flat'], lamPar['gen'], normLL)
    gCFGen.SetName(f'gCFGen{iIter}')
    gCFGen.Write()
    cFit.Write()
    return scattLen, scattLenUnc, normLL, chi2ndf, hCFGen


def ComputeGenCF(args):
    gStyle.SetOptStat(0)

    kStarBW = 50  # MeV/c
    if args.pair == 'DstarK':
        inFileData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaled.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_shifted.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL.root'
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700.root'
        config = '/home/daniel/an/DstarK/cfg_gencf_DstarK_50MeV.yml'

        lightMass = TDatabasePDG.Instance().GetParticle(321).Mass()

        weights1 = [0.78, 0.80, 0.77]
        radii1 = [0.86, 0.95, 0.79]
        radii2 = [2.03, 2.22, 1.91]
    elif args.pair == 'DstarPi':
        inFileData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_mcgp_nopc_kStarBW50MeV_true.root')
        oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}.root'
        if False:
            inFileData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/sbsyst/Distr_mass_2_7_0.20.root')
            oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_sbsyst_mass_2_7_0.20_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}.root'
        if True:
            inFileData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/sbsyst/Distr_mass_2_5_0.07.root')
            oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_sbsyst_mass_2_5_0.07_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}.root'
        config = '/home/daniel/an/DstarPi/cfg_gencf_DstarPi_50MeV.yml'

        lightMass = TDatabasePDG.Instance().GetParticle(211).Mass()

        weights1 = [0.66, 0.69, 0.64]
        radii1 = [0.97, 1.06, 0.89]
        radii2 = [2.52, 2.88, 2.32]
    else:
        print("not implemented")
        sys.exit()

    fitRanges = [[10, 450], [10, 400], [10, 500]]
    fitRanges = [[10, 700], [10, 800], [10, 600]]
    bkgFitRanges = [[300, 1000], [350, 1100], [250, 900]]
    normRange = [1500, 2000]
    heavyMass = TDatabasePDG.Instance().GetParticle(411).Mass()

    # load yaml file with lambda parameter
    with open(config, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    oFile = TFile(oFileName, 'create')
    if oFile.IsZombie():
        sys.exit()

    for comb in ['sc', 'oc']:
        np.random.seed(42)
        if args.pair == 'DstarPi' and comb == 'sc':
            realCoulombFile = TFile('/home/daniel/an/DPi/corrections/Dstar/Dstar_PiplusDplusOutput.root')
        elif args.pair == 'DstarPi' and comb == 'oc':
            realCoulombFile = TFile('/home/daniel/an/DPi/corrections/Dstar/Dstar_PiplusDminusOutput.root')
        elif args.pair == 'DstarK' and comb == 'sc':
            realCoulombFile = TFile('/home/daniel/an/DPi/corrections/Dstar/Dstar_KplusDplusOutput.root')
        elif args.pair == 'DstarK' and comb == 'oc':
            realCoulombFile = TFile('/home/daniel/an/DPi/corrections/Dstar/Dstar_KplusDminusOutput.root')
        fRealCoulomb = realCoulombFile.Get('genuineCF')
        fRealCoulomb.SetLineStyle(9)

        realCoulombFile.Close()

        oFile.mkdir(comb)
        oFile.cd(comb)

        # Compute MC CF
        hSEMC = inFileMC.Get(f'{comb}/SE/sgn/hCharmMassVsKStar0').ProjectionX('hSEMC')
        hMEMC = inFileMC.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hMEMC')

        rebinFactor = round(kStarBW/hSEMC.GetBinWidth(1))
        hSEMC.Rebin(rebinFactor)
        hMEMC.Rebin(rebinFactor)
        hSEMC.Scale(ComputeNormFactor(hSEMC, hMEMC, 1000, 1500))
        hCFMC = hSEMC/hMEMC
        hCFMC.SetName('hCFMC')
        hCFMC.Write()

        regions = fempy.utils.GetRegions(inFileData.Get('sc/SE'))
        nVar = 20 if args.syst else 1
        dCFData = [{} for _ in range(nVar)]
        dSEData = [{} for _ in range(nVar)]
        dhhSEData = [{} for _ in range(nVar)]
        dMEData = [{} for _ in range(nVar)]
        dSEPurity = [{} for _ in range(nVar)]
        dMEPurity = [{} for _ in range(nVar)]
        for iVar in range(nVar):
            for region in regions:
                hhSEData = inFileData.Get(f'{comb}/SE/{region}/hCharmMassVsKStar{iVar}')
                hSEData = hhSEData.ProjectionX(f'hSE_{region}{iVar}')
                dSEData[iVar][region] = hSEData

                hhMEData = inFileData.Get(f'{comb}/ME/{region}/hCharmMassVsKStar{iVar}')
                hMEData = hhMEData.ProjectionX(f'hME_{region}{iVar}')
                dMEData[iVar][region] = hMEData

                rebinFactor = round(kStarBW/hSEData.GetBinWidth(1))
                hSEData.Rebin(rebinFactor)
                hMEData.Rebin(rebinFactor)
                hSEData.Scale(ComputeNormFactor(hSEData, hMEData, normRange[0], normRange[1]))
                hCFData = hSEData/hMEData
                hCFData.SetName(f'hCF_{region}{iVar}')
                dCFData[iVar][region] = hCFData

            dhhSEData[iVar] = inFileData.Get(f'{comb}/SE/hCharmMassVsKStar{iVar}')
            dSEPurity[iVar] = ComputePurity(inFileData.Get(f'{comb}/SE/hCharmMassVsKStar{iVar}'), event='SE', kStarBW=kStarBW) if args.pair == 'DstarPi' else None
            dMEPurity[iVar] = ComputePurity(inFileData.Get(f'{comb}/ME/hCharmMassVsKStar{iVar}'), event='ME', kStarBW=kStarBW) if args.pair == 'DstarPi' else None

        # Compute center of gravity of the bins in the ME
        hGravities = inFileData.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hGravities')
        gGravities = LoadGravities(hGravities, kStarBW)
        gGravities.Write()

        tTrialsStat = TNtuple('tTrialsStat', 'trials', 'a0:a0unc:status:chi2ndf:iVar:purityVar:w1:r1:r2:fitMax:bkgFitMin:bkgFitMax:lFlat:lGen:normLL')

        oFile.mkdir(f'{comb}/stat')
        oFile.cd(f'{comb}/stat')

        hhCFGenStat = TH2D('hhCFGenStat', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)

        if args.pair == 'DstarK':
            lastBin = dhhSEData[0].GetXaxis().FindBin(200*0.9999)
            purity, _ = ComputeIntegratedPurity(dhhSEData[0].ProjectionY(f"Purity_{0}", 1, lastBin))
            lamPar = SumLamPar(LoadLambdaParam(cfg, 0, purity), cfg['treatment'])
        else:
            lamPar = SumLamPar(LoadLambdaParam(cfg), cfg['treatment'])

        for iIter in range(args.bs + 1):  # iter 0 is for the central
            if args.pair == 'DstarK':
                hCFSgn = dCFData[0]['sgn'].Clone(f'hCFSgn{iIter}')
                hCFSbr = dCFData[0]['sbr'].Clone(f'hCFSbr{iIter}')
            elif args.pair == 'DstarPi':
                hSESgn = dSEPurity[0] * dSEData[0]['sgn']
                hMESgn = dMEPurity[0] * dMEData[0]['sgn']
                hCFSgn = hSESgn/hMESgn
            hCFMJ = hCFMC.Clone(f'hCFMC{iIter}')

            scattLen, _, normLL, _, hCFGen = ComputeScattPar(
                sgn=hCFSgn if iIter == 0 else Bootstrap(hCFSgn),
                sbr=None if args.pair == 'DstarPi' else hCFSbr if iIter == 0 else Bootstrap(hCFSbr),
                mj=hCFMJ if iIter == 0 else Bootstrap(hCFMJ),
                iVar=0,
                iIter=iIter,
                comb=comb,
                lamPar=lamPar,
                fitRange=fitRanges[0],
                bkgFitRange=bkgFitRanges[0],
                redMass=RedMass(heavyMass, lightMass) * 1000,
                radius1=radii1[0],
                radius2=radii2[0],
                weights1=weights1[0],
                hhCFGen=hhCFGenStat,
                gGravities=gGravities,
                tTrials=tTrialsStat,
                purityVar=0,
                scalableLL=args.scalableLL,
            )

        oFile.cd(comb)
        tTrialsStat.Write()
        hhCFGenStat.Write()
        hCFGenStat = hCFSgn.Clone('hCFGenStat')
        hCFGenStat.Reset()

        for iBin in range(hCFGenStat.GetNbinsX()):
            hCFGenStatProj = hhCFGenStat.ProjectionY(f'hGenCF_bin{iBin}', iBin+1, iBin+1)
            hCFGenStat.SetBinContent(iBin+1, hCFGenStatProj.GetMean())
            hCFGenStat.SetBinError(iBin+1, hCFGenStatProj.GetStdDev())
        hCFGenStat.Write()

        if args.bs > 0:
            gCFGenStat = ApplyCenterOfGravity(hCFGenStat, gGravities)
        else:
            gCFGenStat = ApplyCenterOfGravity(hCFGen, gGravities)
        
        nLednickyPoints = 500

        if args.syst:
            oFile.mkdir(f'{comb}/tot')
            oFile.cd(f'{comb}/tot')
            hhCFGenTot = TH2D('hhCFGenTot', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)
            tTrialsTot = TNtuple('tTrialsTot', 'trials', 'a0:a0unc:status:chi2ndf:iVar:purityVar:w1:r1:r2:fitMax:bkgFitMin:bkgFitMax:lFlat:lGen:normLL')

            for iIter in range(args.bs):
                iVar = np.random.randint(nVar)

                # non prompt fraction D*: (7.7 +/m 1.3 (stat) +/- 0.2 (syst)) % (syst inherited from pD)
                npFracAbsVar = [0, (0.013**2 + 0.002**2)**0.5, -(0.013**2 + 0.002**2)**0.5][np.random.randint(3)]

                # variation o thermal first
                lightFracRelVar = [0, 0.10, -0.10][np.random.randint(3)]

                # variation o D* meson purity

                if args.pair == 'DstarK':
                    lastBin = dhhSEData[0].GetXaxis().FindBin(200*0.9999)
                    purity, _ = ComputeIntegratedPurity(dhhSEData[iVar].ProjectionY(f"hPurity_{0}", 1, lastBin))
                    DmesonPurityRelVar = [0, 0.02, -0.02][np.random.randint(3)]
                    purity *= (1. + DmesonPurityRelVar)
                    lamPar = SumLamPar(LoadLambdaParam(cfg, npFracAbsVar, purity, lightFracRelVar), cfg['treatment'])

                    if not IsLamParMatValid(LoadLambdaParam(cfg, npFracAbsVar, purity)):
                        fempy.error("lambda parameters don't sum to 1!!")

                    hCFSgn = dCFData[iVar]['sgn'].Clone(f'hCFSgn{iIter}')
                    hCFSbr = dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}')
                    purityVar = 0
                elif args.pair == 'DstarPi':
                    purityVar = np.random.randint(-1, 1)
                    hSESgn = VaryHistogram(dSEPurity[iVar], purityVar) * dSEData[iVar]['sgn']
                    hMESgn = VaryHistogram(dMEPurity[iVar], purityVar) * dMEData[iVar]['sgn']
                    hCFSgn = hSESgn/hMESgn
                    purity = purityVar
                    lamPar = SumLamPar(LoadLambdaParam(cfg, npFracAbsVar, lightFracRelVar=lightFracRelVar), cfg['treatment'])

                radius1, radius2, weight1 = list(zip(radii1, radii2, weights1))[np.random.randint(3)]

                ComputeScattPar(
                    sgn=Bootstrap(hCFSgn.Clone(f'hCFSgn{iIter}')),
                    sbr=None if args.pair == 'DstarPi' else Bootstrap(dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}')),
                    mj=Bootstrap(hCFMC.Clone(f'hCFMC{iIter}')),
                    comb=comb,
                    iVar=iVar,
                    iIter=iIter,
                    lamPar=lamPar,
                    fitRange=fitRanges[np.random.randint(3)],
                    bkgFitRange=bkgFitRanges[np.random.randint(3)],
                    redMass=RedMass(heavyMass, lightMass) * 1000,
                    radius1=radius1,
                    radius2=radius2,
                    weights1=weight1,
                    gGravities=gGravities,
                    hhCFGen=hhCFGenTot,
                    tTrials=tTrialsTot,
                    purityVar=purity,
                    normLL=normLL,
                    scalableLL=args.scalableLL,
                )
            oFile.cd(comb)

            tTrialsTot.Write()
            hhCFGenTot.Write()

            hCFGenTot = hCFSgn.Clone('hCFGenTot')
            for iBin in range(hCFGenTot.GetNbinsX()):
                hCFGenTotProj = hhCFGenTot.ProjectionY(f'hGenCF_bin{iBin}', iBin+1, iBin+1)
                hCFGenTot.SetBinContent(iBin+1, hCFGenTotProj.GetMean())
                hCFGenTot.SetBinError(iBin+1, hCFGenTotProj.GetStdDev())
            hCFGenTot.Write()

            gCFGenSyst = ApplyCenterOfGravity(hCFGenTot, gGravities)
            for iBin in range(60):
                sigma = gCFGenSyst.GetErrorY(iBin)**2 - gCFGenStat.GetErrorY(iBin)**2
                sigma = 0 if sigma < 0 else sigma**0.5
                shift = gCFGenSyst.GetPointY(iBin) - gCFGenStat.GetPointY(iBin)

                gCFGenSyst.SetPoint(iBin, gCFGenSyst.GetPointX(iBin), gCFGenStat.GetPointY(iBin))
                gCFGenSyst.SetPointError(iBin, 0.5*gCFGenSyst.GetErrorX(iBin), (shift**2 + sigma**2)**0.5)
            gCFGenSyst.SetFillColor(38)

            # plot syst curve
            dfTrialsTot = pd.DataFrame(RDataFrame(tTrialsTot).AsNumpy())
            cfVariationTot = [[] for _ in range(nLednickyPoints)]
            for iVar, (scattLen, iVar, r1, r2, w1, normLL) in enumerate(zip(dfTrialsTot['a0'], dfTrialsTot['iVar'], dfTrialsTot['r1'], dfTrialsTot['r2'], dfTrialsTot['w1'], dfTrialsTot['normLL'])):
                fWeightedLLTot = TF1(f"fWeightedLLTot{iVar}", LednickyModel, fitRanges[0][0], fitRanges[0][1], 9)
                fWeightedLLTot.FixParameter(0, r1)
                fWeightedLLTot.FixParameter(1, r2)
                fWeightedLLTot.FixParameter(2, w1)
                fWeightedLLTot.FixParameter(3, scattLen)
                fWeightedLLTot.FixParameter(4, 0.)
                fWeightedLLTot.FixParameter(5, 0.)
                fWeightedLLTot.FixParameter(6, RedMass(lightMass, heavyMass)*1000)
                fWeightedLLTot.FixParameter(7, 1 if comb == 'sc' else -1)
                fWeightedLLTot.FixParameter(8, normLL)
                for iPoint in range(nLednickyPoints):
                    cfVariationTot[iPoint].append(fWeightedLLTot.Eval(float(iPoint)/nLednickyPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0]))

        cFinalFit = TCanvas(f'cFinalFit_{comb}', '', 600, 600)
        cFinalFit.DrawFrame(0, 0, 500, 2, fempy.utils.format.TranslateToLatex(';__kStarMeV__;__C__'))
        gCFGenStat.SetMarkerStyle(33)
        gCFGenStat.SetMarkerSize(2)

        # plot lednicky curves
        dfTrials = pd.DataFrame(RDataFrame(tTrialsStat).AsNumpy())
        cfVariationStat = [[] for _ in range(nLednickyPoints)]
        for iIter, (scattLen, normLL) in enumerate(zip(dfTrials['a0'], dfTrials['normLL'])):
            fWeightedLLStat = TF1(f"fWeightedLLStat{iIter}", LednickyModel, fitRanges[0][0], fitRanges[0][1], 9)
            fWeightedLLStat.FixParameter(0, radii1[0])
            fWeightedLLStat.FixParameter(1, radii2[0])
            fWeightedLLStat.FixParameter(2, weights1[0])
            fWeightedLLStat.FixParameter(3, scattLen)
            fWeightedLLStat.FixParameter(4, 0.)
            fWeightedLLStat.FixParameter(5, 0.)
            fWeightedLLStat.FixParameter(6, RedMass(lightMass, heavyMass)*1000)
            fWeightedLLStat.FixParameter(7, 1 if comb == 'sc' else -1)
            fWeightedLLStat.FixParameter(8, normLL)

            for iPoint in range(nLednickyPoints):
                cfVariationStat[iPoint].append(fWeightedLLStat.Eval(float(iPoint)/nLednickyPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0]))

        # Compute Scat param
        hScatParStat = TH1D('hScatParStat', ';a_{0} (fm);Counts', 200, -1, 1)
        tTrialsStat.Project('hScatParStat', 'a0')
        scatPar = hScatParStat.GetMean()
        scatParStatUnc = hScatParStat.GetStdDev()
        if args.syst:
            gLLTot = TGraphErrors(1)
            for iPoint in range(nLednickyPoints):
                gLLTot.SetPoint(iPoint, float(iPoint)/nLednickyPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0], np.average(cfVariationStat[iPoint]))
                shift = np.average(cfVariationStat[iPoint]) - np.average(cfVariationTot[iPoint])
                sigma = np.std(cfVariationTot[iPoint])
                gLLTot.SetPointError(iPoint, 0, (shift**2 + sigma**2)**0.5)
            gLLTot.SetFillColor(42)

            hScatParTot = TH1D('hScatParTot', ';a_{0} (fm);Counts', 200, -1, 1)
            tTrialsTot.Project('hScatParTot', 'a0')
            scatParSystUnc = (hScatParTot.GetStdDev()**2 - scatParStatUnc**2)**0.5

        gLLStat = TGraphErrors(1)
        for iPoint in range(nLednickyPoints):
            gLLStat.SetPoint(iPoint, float(iPoint)/nLednickyPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0], np.average(cfVariationStat[iPoint]))
            gLLStat.SetPointError(iPoint, 0, np.std(cfVariationStat[iPoint]))
        gLLStat.SetFillColor(46)
        gLLStat.SetLineColor(46)

        if args.bs > 0 and args.syst:
            chi2 = 0
            # Chi2 recalculation
            nPoints = 9
            ndf = nPoints - 1
            for iPoint in range(ndf):
                kStar = gCFGenStat.GetPointX(iPoint)
                cfFit = gLLStat.Eval(kStar)

                cfData = gCFGenStat.GetPointY(iPoint)
                cfDataStatUnc = gCFGenStat.GetErrorY(iPoint)
                cfDataSystUnc = gCFGenSyst.GetErrorY(iPoint) if args.syst else 0

                cfDataTotUnc = (cfDataStatUnc**2 + cfDataSystUnc**2)**0.5
                chi2 += ((cfData - cfFit)/cfDataTotUnc)**2

        step = 0.05
        tl = TLatex()
        tl.SetNDC()
        tl.SetTextSize(0.035)
        tl.SetTextFont(42)
        tl.DrawLatex(0.2, 0.85, fempy.utils.format.TranslateToLatex(f'k{args.pair}_{comb}'))
        tl.DrawLatex(0.2, 0.85 - step, f'a_{{0}} = {scatPar:.3f} #pm {scatParStatUnc:.3f} (stat){f" #pm {scatParSystUnc:.3f} (syst)" if args.syst else ""}')
        if args.bs > 0 and args.syst:
            tl.DrawLatex(0.2, 0.85 - 2 * step, f'#chi^{{2}}/ndf = {chi2:.0f} / {ndf:.0f}')

        fCoulombLL = TF1("fCoulombLL", LednickyModel, fitRanges[0][0], fitRanges[0][1], 9)
        fCoulombLL.FixParameter(0, radii1[0])
        fCoulombLL.FixParameter(1, radii2[0])
        fCoulombLL.FixParameter(2, weights1[0])
        fCoulombLL.FixParameter(3, 0)
        fCoulombLL.FixParameter(4, 0.)
        fCoulombLL.FixParameter(5, 0.)
        fCoulombLL.FixParameter(6, RedMass(lightMass, heavyMass)*1000)
        fCoulombLL.FixParameter(7, 1 if comb == 'sc' else -1)
        fCoulombLL.FixParameter(8, 1)
        fCoulombLL.SetLineColor(kBlue)

        leg = TLegend(0.6, 0.75, .9, 0.9)
        if args.syst:
            leg.AddEntry(gCFGenSyst, 'Data', 'pef')
        else:
            leg.AddEntry(gCFGenStat, 'Data', 'pef')
        leg.AddEntry(gLLStat, 'Fit LL', 'f')
        leg.AddEntry(fCoulombLL, 'LL Coulomb', 'l')
        leg.AddEntry(fRealCoulomb, 'Real Coulomb', 'l')

        leg.Draw()
        if args.syst:
            gLLTot.Draw('same e3')
        gLLStat.Draw('same e3')
        fCoulombLL.Draw('same')
        fRealCoulomb.Draw('same')
        if args.syst:
            gCFGenSyst.Draw('same pe2')
        else:
            gCFGenStat.Draw('same pe')

        gCFGenStat.Draw('same pe')

        cFinalFit.Modified()
        cFinalFit.Update()
        cFinalFit.Write()
        gCFGenStat.SetName('gCFGenStat')
        gCFGenStat.Write()
        if args.syst:
            gCFGenSyst.SetName('gCFGenSyst')
            gCFGenSyst.Write()
            gFit = gLLTot.Clone()
            gFit.SetName('gFit')
            gFit.Write()

        gBrackets = ComputeBinBrackets(hCFGenStat)
        gBrackets.Write()

    print(f'output saved in {oFileName}')
    oFile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarPi', 'DstarK'))
    parser.add_argument('--syst', action='store_true', default=False)
    parser.add_argument('--scalableLL', action='store_true', default=False)
    parser.add_argument('--shift', action='store_true', default=False)
    parser.add_argument('--bs', type=int, default=0)
    args = parser.parse_args()

    LednickyModel = ScalableShifted300GeneralCoulombLednickyTwoRadii if args.shift else ScalableGeneralCoulombLednickyTwoRadii

    ComputeGenCF(args)
