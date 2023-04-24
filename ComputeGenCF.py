import sys
import yaml
import random
import argparse

import pandas as pd
import numpy as np

from ROOT import TFile, TF1, TGraphErrors, TDatabasePDG, gInterpreter, TCanvas, kBlue, TLatex, gStyle, TLegend, gRandom, TNtuple, TH2D, RDataFrame, TH1D
gInterpreter.ProcessLine('#include "fempy/utils/functions.h"')
gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter, GeneralCoulombLednicky

import fempy


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


def LoadLambdaParam(cfgCentr, npFracVar=1.):
    cfgVar = dict(cfgCentr)
    cfgVar['heavy'][1]['nonprompt']['frac'] *= npFracVar
    cfgVar['heavy'][0]['prompt']['frac'] = 1 - cfgVar['heavy'][1]['nonprompt']['frac']
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


def ComputeScattPar(**kwargs):
    hCFSgn = kwargs['sgn']
    hCFSbr = kwargs['sbr']
    hCFMJ = kwargs['mj']
    lamPar = kwargs['lamPar']
    gGravities = kwargs['gGravities']
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
    fCoulomb = TF1("fCoulomb", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
    fCoulomb.FixParameter(0, radius1)
    fCoulomb.FixParameter(1, radius2)
    fCoulomb.FixParameter(2, weight1)
    fCoulomb.FixParameter(3, 0)
    fCoulomb.FixParameter(4, 0.)
    fCoulomb.FixParameter(5, 0.)
    fCoulomb.FixParameter(6, 1 if comb == 'sc' else -1)
    fCoulomb.FixParameter(7, redMass)
    fCoulomb.SetLineColor(kBlue)

    # Fit the CF
    gCFGen = ApplyCenterOfGravity(hCFGen, gGravities)
    fWeightedLL = TF1(f"fWeightedLL{iIter}", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
    fWeightedLL.FixParameter(0, radius1)
    fWeightedLL.FixParameter(1, radius2)
    fWeightedLL.FixParameter(2, weight1)
    fWeightedLL.SetParameter(3, 0.1)
    fWeightedLL.SetParLimits(3, -1, 1)
    fWeightedLL.FixParameter(4, 0.)
    fWeightedLL.FixParameter(5, 0.)
    fWeightedLL.FixParameter(6, 1 if comb == 'sc' else -1)
    fWeightedLL.FixParameter(7, redMass)

    status = gCFGen.Fit(fWeightedLL, "SMRQ+0").Status()
    chi2ndf = fWeightedLL.GetChisquare()/fWeightedLL.GetNDF()
    scattLen = fWeightedLL.GetParameter(3)
    scattLenUnc = fWeightedLL.GetParError(3)

    # Draw canvas
    cFit = TCanvas(f'cFit_{comb}{iIter}', '', 600, 600)
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
        tTrials.Fill(scattLen, scattLenUnc, status, chi2ndf, iVar, purityVar, weight1, radius1, radius2, fitRange[1], bkgFitRange[0], bkgFitRange[1], lamPar['flat'], lamPar['gen'])
    cFit.Write()


def ComputeGenCF(args):
    random.seed(42)
    gStyle.SetOptStat(0)

    kStarBW = 50  # MeV/c
    if args.pair == 'DstarK':
        inFileData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}.root'
        config = '/home/daniel/an/DstarK/cfg_gencf_DstarK_50MeV.yml'

        lightMass = TDatabasePDG.Instance().GetParticle(321).Mass()

        weights1 = [0.78, 0.80, 0.77]
        radii1 = [0.86, 0.95, 0.79]
        radii2 = [2.03, 2.22, 1.91]
    elif args.pair == 'DstarPi':
        fitRanges = [[10, 800], [10, 600], [10, 1000]]
        inFileData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_mcgp_nopc_kStarBW50MeV_true.root')
        oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}.root'
        config = '/home/daniel/an/DstarPi/cfg_gencf_DstarPi_50MeV.yml'

        lightMass = TDatabasePDG.Instance().GetParticle(211).Mass()
        
        weights1 = [0.66, 0.69, 0.64]
        radii1 = [0.97, 1.06, 0.89]
        radii2 = [2.52, 2.88, 2.32]
    else:
        print("not implemented")
        sys.exit()

    fitRanges = [[10, 450], [10, 400], [10, 500]]
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

    # Variations on lambda parameters
    lamParCentr = SumLamPar(LoadLambdaParam(cfg), cfg['treatment'])
    lamParSharp = SumLamPar(LoadLambdaParam(cfg, 1.1), cfg['treatment'])
    lamParFlat = SumLamPar(LoadLambdaParam(cfg, 0.9), cfg['treatment'])

    lamPar = lamParCentr

    oFile = TFile(oFileName, 'recreate')
    for comb in ['sc', 'oc']:
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

            dSEPurity[iVar] = ComputePurity(inFileData.Get(f'{comb}/SE/hCharmMassVsKStar{iVar}'), event='SE', kStarBW=kStarBW) if args.pair == 'DstarPi' else None
            dMEPurity[iVar] = ComputePurity(inFileData.Get(f'{comb}/ME/hCharmMassVsKStar{iVar}'), event='ME', kStarBW=kStarBW) if args.pair == 'DstarPi' else None

        # Compute center of gravity of the bins in the ME
        hGravities = inFileData.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hGravities')
        gGravities = LoadGravities(hGravities, kStarBW)
        gGravities.Write()

        tTrialsStat = TNtuple('tTrialsStat', 'trials', 'a0:a0unc:status:chi2ndf:iVar:purityVar:w1:r1:r2:fitMax:bkgFitMin:bkgFitMax:lFlat:lGen')

        oFile.mkdir(f'{comb}/stat')
        oFile.cd(f'{comb}/stat')

        hhCFGenStat = TH2D('hhCFGenStat', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)

        for iIter in range(args.bs + 1):  # iter 0 is for the central
            if args.pair == 'DstarK':
                hCFSgn = dCFData[0]['sgn'].Clone(f'hCFSgn{iIter}')
                hCFSbr = dCFData[0]['sbr'].Clone(f'hCFSbr{iIter}')
            elif args.pair == 'DstarPi':
                hSESgn = dSEPurity[0] * dSEData[0]['sgn']
                hMESgn = dMEPurity[0] * dMEData[0]['sgn']
                hCFSgn = hSESgn/hMESgn
            hCFMJ = hCFMC.Clone(f'hCFMC{iIter}')

            ComputeScattPar(
                sgn=hCFSgn if iIter == 0 else Bootstrap(hCFSgn),
                sbr=None if args.pair=='DstarPi' else hCFSbr if iIter == 0 else Bootstrap(hCFSbr),
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

        gCFGenStat = ApplyCenterOfGravity(hCFGenStat, gGravities)
        cFinalFit = TCanvas(f'cFinalFit_{comb}', '', 600, 600)
        cFinalFit.DrawFrame(0, 0, 500, 2, fempy.utils.format.TranslateToLatex(';__kStar_MeV__;__C__'))
        leg = TLegend(0.6, 0.75, .9, 0.9)

        if args.syst:
            oFile.mkdir(f'{comb}/tot')
            oFile.cd(f'{comb}/tot')
            hhCFGenTot = TH2D('hhCFGenTot', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)
            tTrialsTot = TNtuple('tTrialsTot', 'trials', 'a0:a0unc:status:chi2ndf:iVar:purityVar:w1:r1:r2:fitMax:bkgFitMin:bkgFitMax:lFlat:lGen')

            for iIter in range(args.bs):
                iVar = random.choice(range(0, nVar))
                if args.pair == 'DstarK':
                    hCFSgn = dCFData[iVar]['sgn'].Clone(f'hCFSgn{iIter}')
                    hCFSbr = dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}')
                elif args.pair == 'DstarPi':
                    purityVar = random.choice([0, +1, -1])
                    hSESgn = VaryHistogram(dSEPurity[iVar], purityVar) * dSEData[iVar]['sgn']
                    hMESgn = VaryHistogram(dMEPurity[iVar], purityVar) * dMEData[iVar]['sgn']
                    hCFSgn = hSESgn/hMESgn
                radius1, radius2, weight1 = random.choice(list(zip(radii1, radii2, weights1)))

                ComputeScattPar(
                    sgn=Bootstrap(hCFSgn.Clone(f'hCFSgn{iIter}')),
                    sbr=None if args.pair=='DstarPi' else Bootstrap(dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}')),
                    mj=Bootstrap(hCFMC.Clone(f'hCFMC{iIter}')),
                    comb=comb,
                    iVar=iVar,
                    iIter=iIter,
                    lamPar=random.choice([lamParCentr, lamParFlat, lamParSharp]),
                    fitRange=random.choice(fitRanges),
                    bkgFitRange=random.choice(bkgFitRanges),
                    redMass=RedMass(heavyMass, lightMass) * 1000,
                    radius1=radius1,
                    radius2=radius2,
                    weights1=weight1,
                    gGravities=gGravities,
                    hhCFGen=hhCFGenTot,
                    tTrials=tTrialsTot,
                    purityVar=purityVar,
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

            gCFGenTot = ApplyCenterOfGravity(hCFGenTot, gGravities)
            for iBin in range(60):
                syst = gCFGenTot.GetErrorY(iBin)**2 - gCFGenStat.GetErrorY(iBin)**2
                syst = 0 if syst < 0 else syst**0.5
                gCFGenTot.SetPointError(iBin, 0.5*gCFGenTot.GetErrorX(iBin), syst)

            gCFGenTot.SetFillColor(38)
            gCFGenTot.Draw('same pe2')
            gCFGenStat.Draw('same pe')
        
            leg.AddEntry(gCFGenTot, 'Data')
            leg.Draw()
            
            # plot lednicky curves
            nPoints = 500

            dfTrials = pd.DataFrame(RDataFrame(tTrialsStat).AsNumpy())
            cfVariationStat = [[] for _ in range(nPoints)]
            for iIter, scattLen in enumerate(dfTrials['a0']):
                fWeightedLLStat = TF1(f"fWeightedLLStat{iIter}", WeightedCoulombLednicky, fitRanges[0][0], fitRanges[0][1], 8)
                fWeightedLLStat.FixParameter(0, radii1[0])
                fWeightedLLStat.FixParameter(1, radii2[0])
                fWeightedLLStat.FixParameter(2, weights1[0])
                fWeightedLLStat.FixParameter(3, scattLen)
                fWeightedLLStat.FixParameter(4, 0.)
                fWeightedLLStat.FixParameter(5, 0.)
                fWeightedLLStat.FixParameter(6, 1 if comb == 'sc' else -1)
                fWeightedLLStat.FixParameter(7, RedMass(lightMass, heavyMass)*1000)
                
                for iPoint in range(nPoints):
                    cfVariationStat[iPoint].append(fWeightedLLStat.Eval(float(iPoint)/nPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0]))

            gLLStat = TGraphErrors(1)
            for iPoint in range(nPoints):
                gLLStat.SetPoint(iPoint, float(iPoint)/nPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0], np.average(cfVariationStat[iPoint]))
                gLLStat.SetPointError(iPoint, 0, np.std(cfVariationStat[iPoint]))
            gLLStat.SetFillColor(46)

            # plot syst curve
            dfTrialsTot = pd.DataFrame(RDataFrame(tTrialsTot).AsNumpy())
            cfVariationTot = [[] for _ in range(nPoints)]
            for iVar, (scattLen, iVar, r1, r2, w1) in enumerate(zip(dfTrialsTot['a0'], dfTrialsTot['iVar'], dfTrialsTot['r1'], dfTrialsTot['r2'], dfTrialsTot['w1'])):
                fWeightedLLTot = TF1(f"fWeightedLLTot{iVar}", WeightedCoulombLednicky, fitRanges[0][0], fitRanges[0][1], 8)
                fWeightedLLTot.FixParameter(0, r1)
                fWeightedLLTot.FixParameter(1, r2)
                fWeightedLLTot.FixParameter(2, w1)
                fWeightedLLTot.FixParameter(3, scattLen)
                fWeightedLLTot.FixParameter(4, 0.)
                fWeightedLLTot.FixParameter(5, 0.)
                fWeightedLLTot.FixParameter(6, 1 if comb == 'sc' else -1)
                fWeightedLLTot.FixParameter(7, RedMass(lightMass, heavyMass)*1000)
                
                for iPoint in range(nPoints):
                    cfVariationTot[iPoint].append(fWeightedLLTot.Eval(float(iPoint)/nPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0]))

            gLLSyst = TGraphErrors(1)
            for iPoint in range(nPoints):
                gLLSyst.SetPoint(iPoint, float(iPoint)/nPoints*(fitRanges[0][1] - fitRanges[0][0]) + fitRanges[0][0], np.average(cfVariationTot[iPoint]))
                gLLSyst.SetPointError(iPoint, 0, np.std(cfVariationTot[iPoint]))
            gLLSyst.SetFillColor(42)
            gLLSyst.Draw('same e3')
            gLLStat.Draw('same e3')

            # Compute Scat param
            hScatParStat = TH1D('hScatParStat', ';a_{0} (fm);Counts', 200, -1, 1)
            tTrialsStat.Project('hScatParStat', 'a0')
            scatPar = hScatParStat.GetMean()
            scatParStatUnc = hScatParStat.GetStdDev()
            
            hScatParTot = TH1D('hScatParTot', ';a_{0} (fm);Counts', 200, -1, 1)
            tTrialsTot.Project('hScatParTot', 'a0')
            scatParSystUnc = (hScatParTot.GetStdDev()**2 - scatParStatUnc**2)**0.5
            
            step = 0.05
            tl = TLatex()
            tl.SetNDC()
            tl.SetTextSize(0.035)
            tl.SetTextFont(42)
            tl.DrawLatex(0.2, 0.85, fempy.utils.format.TranslateToLatex(f'k{args.pair}_{comb}'))
            tl.DrawLatex(0.2, 0.85 - step, f'a_{{0}} = {scatPar:.3f} #pm {scatParStatUnc:.3f} (stat) #pm {scatParSystUnc:.3f} (syst)')

            cFinalFit.Modified()
            cFinalFit.Update()
            cFinalFit.Write()
    print(f'output saved in {oFileName}')
    oFile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarPi', 'DstarK'))
    parser.add_argument('--syst', action='store_true', default=False)
    parser.add_argument('--bs', type=int, default=0)
    args = parser.parse_args()

    ComputeGenCF(args)
