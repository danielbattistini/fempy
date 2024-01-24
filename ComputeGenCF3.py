import sys
import yaml
import argparse
import copy

import pandas as pd
import numpy as np

from ROOT import TFile, TF1, TGraphErrors, gInterpreter, gStyle, gRandom, TH2D, TH1D
gInterpreter.ProcessLine('#include "combfit/functions.h"')
gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter

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
    # cPurity = TCanvas('cPurity', '', 600, 600)
    # cPurity.SaveAs(f'cPurity_{hist.GetName()}_{event}.pdf[')
    for iBin in range(nBins):
        firstBin = hist.FindBin(iBin * kStarBW * 1.0001)
        lastBin = hist.FindBin((iBin + 1) * kStarBW * 0.9999)

        hCharmMass = hist.ProjectionY('', firstBin, lastBin)
        hCharmMass.SetTitle(f'k* bin = {iBin}')
        fitter = MassFitter(hCharmMass, 'gaus', 'powex', 0.141, 0.154)
        fitter.SetFitSettings(f'413_gaus_powex_DstarFemto_{event}')
        fitter.Fit()
        # fitter.Draw(cPurity, 'data_minus_bkg')
        # cPurity.SaveAs(f'cPurity_{hist.GetName()}_{event}.pdf')

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
    # cPurity.SaveAs(f'cPurity_{hist.GetName()}_{event}.pdf]')

    return hPurity


def VaryHistogram(hist, method):
    if method == 'bs':
        return Bootstrap(hist)
    elif isinstance(method, int):
        hVaried = hist.Clone()
        hVaried.Reset()
        for iBin in range(hVaried.GetNbinsX()):
            hVaried.SetBinContent(iBin+1, hist.GetBinContent(iBin+1) + method * hist.GetBinError(iBin+1))
            hVaried.SetBinError(iBin+1, hist.GetBinError(iBin+1))
        return hVaried
    else:
        fempy.error("not implemented")

from ROOT import TCanvas

def ComputeIntegratedPurity(hMass, fitRange=[0.141, 0.154], name=''):
    cPurity = TCanvas('cPurity', '', 600, 600)

    firstBin = hMass.GetXaxis().FindBin(fitRange[0]*1.0001)
    lastBin = hMass.GetXaxis().FindBin(fitRange[1]*0.9999)

    particleYield = hMass.Integral(firstBin, lastBin)

    fitter = MassFitter(hMass, 'gaus', 'powex', fitRange[0], fitRange[1])
    fitter.SetFitSettings(f'413_gaus_powex_DstarKFemto_SE')
    fitter.Fit()
    fitter.Draw(cPurity, 'data_minus_bkg')
    cPurity.SaveAs(f'cPurity_{hMass.GetName()}_{name}.pdf')

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
    hhCFGen = kwargs['hhCFGen']
    iIter = kwargs['iIter']
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
        # hCFNorm = (hCFSgn - lamPar['sb'] * hCFSbr) / (hCFMJ * (lamPar['gen'] + lamPar['flat']))
        # hCFNorm.SetName(f'hCFNorm{iIter}')
        # hCFNorm.Write()

        # # Fit the baseline
        # fBaseLine = TF1(f'fBaseLine{iIter}', '[0]', 0, 3000)
        # ApplyCenterOfGravity(hCFNorm, gGravities).Fit(fBaseLine, 'Q', '', bkgFitRange[0], bkgFitRange[1])
        # blNorm = fBaseLine.GetParameter(0)

        # Compute the total background mode
        def bkgModel(x, par):
            bin = hCFSbr.FindBin(x[0])

            sbr = hCFSbr.GetBinContent(bin)
            c_non_femto = par[0] * (hCFMJ.GetBinContent(bin) + par[1] * x[0])
            return lamPar['sb'] * sbr + c_non_femto * (lamPar['gen'] + lamPar['flat'])

        fCFBkgModel = TF1('fCFBkgModel', bkgModel, bkgFitRange[0], bkgFitRange[1], 2) 
        fCFBkgModel.SetParLimits(0, 0, 2)
        fCFBkgModel.SetParLimits(1, -0.1, 0.1)
        ApplyCenterOfGravity(hCFSgn, gGravities).Fit(fCFBkgModel, 'MR+', '', bkgFitRange[0], bkgFitRange[1])
        # hCFBkg = lamPar['sb'] * hCFSbr + blNorm * hCFMJ * (lamPar['gen'] + lamPar['flat'])
        # hCFBkg.SetName(f'hCFBkg{iIter}')
        # hCFBkg.Write()

        # Compute the Gen CF
        
        hCFGen = hCFSgn.Clone(f'hCFGen{iIter}')
        hCFGen.Reset()
        # hCFFlat = MakeFlatHist(hCFSgn, lamPar['flat'], 'hCFFlat')
        for iBin in range(hCFGen.GetNbinsX()):
            kStar = hCFGen.GetBinCenter(iBin + 1)
            sgn = hCFSgn.GetBinContent(iBin + 1)
            sbr = hCFSbr.GetBinContent(iBin + 1)
            mj = hCFMJ.GetBinContent(iBin + 1)
            non_femto = (fCFBkgModel.GetParameter(0) * (mj + kStar * fCFBkgModel.GetParameter(1)))
            hCFGen.SetBinContent(iBin + 1, ((sgn - lamPar['sb'] * sbr)/non_femto - lamPar['flat'])/lamPar['gen'])
            hCFGen.SetBinError(iBin + 1, 0) # Error determined with the bootstrap method

        # hCFGen = (hCFSgn - lamPar['sb'] * hCFSbr)/(blNorm * hCFMJ) - hCFFlat
        # hCFGen.Scale(1./lamPar['gen'])
        # hCFGen.SetName(f'hCFGen{iIter}')
        # hCFGen.Write()

    # Add the CF to the hist2D
    if iIter > 0:
        for iBin in range(hCFGen.GetNbinsX()):
            hhCFGen.Fill(hCFGen.GetBinCenter(iBin+1), hCFGen.GetBinContent(iBin+1))

    gCFGen = ApplyCenterOfGravity(hCFGen, gGravities)
    gCFGen.SetName(f'gCFGen{iIter}')
    gCFGen.Write()
    return hCFGen

normRange = [1500, 2000]

def ComputeGenCF(args):
    gStyle.SetOptStat(0)

    kStarBW = 50  # MeV/c
    if args.pair == 'DstarK':
        inFileData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFDebug_nopc_kStarBW50MeV_fromq_bs{args.bs}{"syst" if args.syst else ""}.root'
        config = '/home/daniel/an/DstarK/cfg_gencf_DstarK_50MeV.yml'
        bkgFitRanges = [[250, 500], [200, 400], [300, 600]]
    elif args.pair == 'DstarPi':
        inFileData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_mcgp_nopc_kStarBW50MeV_true.root')
        oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFDebug_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}.root'
        oFileName = f'/home/daniel/an/DstarPi/20_luuksel/GenCFDebug_nopc_kStarBW50MeV_bs{args.bs}{"syst" if args.syst else ""}_bkgFitRange300-1000.root'
        config = '/home/daniel/an/DstarPi/cfg_gencf_DstarPi_50MeV.yml'
        bkgFitRanges = [[300, 1000], [250, 900], [350, 1100]]
    else:
        print("not implemented")
        sys.exit()

    oFile = TFile(oFileName, 'create')
    if oFile.IsZombie():
        sys.exit()

    # load yaml file with lambda parameter
    with open(config, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    for comb in ['sc', 'oc']:
        np.random.seed(42)

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

        oFile.mkdir(f'{comb}/stat')
        oFile.cd(f'{comb}/stat')

        hhCFGenStat = TH2D('hhCFGenStat', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)

        if args.pair == 'DstarK':
            lastBin = dhhSEData[0].GetXaxis().FindBin(200*0.9999)
            purity, _ = ComputeIntegratedPurity(dhhSEData[0].ProjectionY(f"Purity_{0}", 1, lastBin), name=f'{comb}_SE')
            lamPar = SumLamPar(LoadLambdaParam(cfg, 0, purity), cfg['treatment'])
        else:
            lamPar = SumLamPar(LoadLambdaParam(cfg), cfg['treatment'])

        for iIter in range(args.bs):
            if args.pair == 'DstarK':
                hCFSgn = dCFData[0]['sgn'].Clone(f'hCFSgn{iIter}')
                hCFSbr = dCFData[0]['sbr'].Clone(f'hCFSbr{iIter}')
            elif args.pair == 'DstarPi':
                hSESgn = dSEPurity[0] * dSEData[0]['sgn']
                hMESgn = dMEPurity[0] * dMEData[0]['sgn']
                hCFSgn = hSESgn/hMESgn
            hCFMJ = hCFMC.Clone(f'hCFMC{iIter}')

            ComputeScattPar(
                sgn=Bootstrap(hCFSgn),
                sbr=None if args.pair == 'DstarPi' else Bootstrap(hCFSbr),
                mj=Bootstrap(hCFMJ),
                iVar=0,
                iIter=iIter,
                lamPar=lamPar,
                bkgFitRange=bkgFitRanges[0],
                hhCFGen=hhCFGenStat,
                gGravities=gGravities,
            )

        oFile.cd(comb)
        hhCFGenStat.Write()
        
        hCFGenStat = hCFSgn.Clone('hCFGenStat')
        hCFGenStat.Reset()
        for iBin in range(hCFGenStat.GetNbinsX()):
            hCFGenStatProj = hhCFGenStat.ProjectionY(f'hGenCF_bin{iBin}', iBin+1, iBin+1)
            hCFGenStat.SetBinContent(iBin+1, hCFGenStatProj.GetMean())
            hCFGenStat.SetBinError(iBin+1, hCFGenStatProj.GetStdDev())
        hCFGenStat.Write()

        gCFGenStat = ApplyCenterOfGravity(hCFGenStat, gGravities)
        gCFGenStat.SetName('gCFGenStat')
        gCFGenStat.Write()

        if args.syst:
            oFile.mkdir(f'{comb}/tot')
            oFile.cd(f'{comb}/tot')

            hhCFGenTot = TH2D('hhCFGenTot', '', dCFData[0]['sgn'].GetNbinsX(), 0, 3000, 2000, 0, 2)


            for iIter in range(args.bs):
                iVar = np.random.randint(nVar)

                # non prompt fraction D*: (7.7 +/m 1.3 (stat) +/- 0.2 (syst)) % (syst inherited from pD)
                npFracAbsVar = [0, (0.013**2 + 0.002**2)**0.5, -(0.013**2 + 0.002**2)**0.5][np.random.randint(3)]
                # npFracAbsVar = 0 #!

                # variation o thermal first
                lightFracRelVar = [0, 0.10, -0.10][np.random.randint(3)]
                # lightFracRelVar = 0 #!
                
                hCFMJ = hCFMC.Clone(f'hCFMC{iIter}')
                if args.pair == 'DstarK': #!
                    lastBin = dhhSEData[0].GetXaxis().FindBin(200*0.9999)
                    purity, _ = ComputeIntegratedPurity(dhhSEData[iVar].ProjectionY(f"hPurity_{iVar}", 1, lastBin), name=f'{comb}_{iIter}_{iVar}')
                    DmesonPurityRelVar = [0, 0.02, -0.02][np.random.randint(3)]
                    # DmesonPurityRelVar = 0 #!
                    purity *= (1. + DmesonPurityRelVar)
                    lamPar = SumLamPar(LoadLambdaParam(cfg, npFracAbsVar, purity, lightFracRelVar), cfg['treatment'])


                    
                    
                    
                    # if not IsLamParMatValid(LoadLambdaParam(cfg, npFracAbsVar, purity)):
                    #     fempy.error("lambda parameters don't sum to 1!!")

                    # hCFSgn = Bootstrap(dCFData[iVar]['sgn'].Clone(f'hCFSgn{iIter}'))
                    # hCFSbr = Bootstrap(dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}'))

                    # lastBin = dhhSEData[0].GetXaxis().FindBin(200*0.9999)
                    # purity, _ = ComputeIntegratedPurity(dhhSEData[0].ProjectionY(f"Purity_{0}", 1, lastBin), name=f'{comb}_SE')
                    # lamPar = SumLamPar(LoadLambdaParam(cfg, 0, purity), cfg['treatment'])

                    hCFSgn = dCFData[iVar]['sgn'].Clone(f'hCFSgn{iIter}')
                    hCFSbr = dCFData[iVar]['sbr'].Clone(f'hCFSbr{iIter}')
                    ComputeScattPar(
                        sgn=Bootstrap(hCFSgn),
                        sbr=Bootstrap(hCFSbr),
                        mj=Bootstrap(hCFMJ),
                        iVar=0,
                        iIter=iIter,
                        lamPar=lamPar,
                        bkgFitRange=bkgFitRanges[np.random.randint(3)], #!
                        hhCFGen=hhCFGenTot,
                        gGravities=gGravities,
                    )
                elif args.pair == 'DstarPi':
                    lamPar = SumLamPar(LoadLambdaParam(cfg, npFracAbsVar, lightFracRelVar=lightFracRelVar), cfg['treatment'])
                    hSESgn = dSEPurity[iVar] * dSEData[iVar]['sgn']
                    hMESgn = dMEPurity[iVar] * dMEData[iVar]['sgn']

                    # purityVar = np.random.randint(-1, 2)
                    # purityVar = np.random.randint(-1, 2)

                    hSESgn = VaryHistogram(dSEPurity[iVar], np.random.randint(-1, 2)) * dSEData[iVar]['sgn']
                    hMESgn = VaryHistogram(dMEPurity[iVar], np.random.randint(-1, 2)) * dMEData[iVar]['sgn']
                    hSESgn.Scale(ComputeNormFactor(hSESgn, hMESgn, normRange[0], normRange[1]))
                    hCFSgn = hSESgn/hMESgn

                    ComputeScattPar(
                        sgn=Bootstrap(hCFSgn), #!
                        sbr=None,
                        mj=Bootstrap(hCFMJ),
                        iVar=iVar,
                        iIter=iIter,
                        lamPar=lamPar,
                        bkgFitRange=bkgFitRanges[np.random.randint(3)],
                        hhCFGen=hhCFGenTot,
                        gGravities=gGravities,
                    )

            oFile.cd(comb)
            hhCFGenTot.Write()
            
            hCFGenTot = hCFSgn.Clone('hCFGenTot')
            hCFGenTot.Reset()
            hCFGenSyst = hCFSgn.Clone('hCFGenSyst')
            hCFGenSyst.Reset()
            for iBin in range(hCFGenTot.GetNbinsX()):
                hCFGenTotProj = hhCFGenTot.ProjectionY(f'hGenCF_bin{iBin}', iBin+1, iBin+1)

                muTot = hCFGenTotProj.GetMean()
                sigmaTot = hCFGenTotProj.GetStdDev()

                # to remove outliers
                hCFGenTotProj.GetXaxis().SetRangeUser(muTot - 4 * sigmaTot, muTot + 4 * sigmaTot)
                
                # recalculation
                muTot = hCFGenTotProj.GetMean()
                sigmaTot = hCFGenTotProj.GetStdDev()

                muStat = hCFGenStat.GetBinContent(iBin + 1)
                sigmaStat = hCFGenStat.GetBinError(iBin + 1)
                shift = muStat - muTot

                hCFGenTot.SetBinContent(iBin+1, muTot)
                hCFGenTot.SetBinError(iBin+1, sigmaTot)

                sigmaSyst = sigmaTot ** 2 - sigmaStat ** 2 + shift ** 2
                sigmaSyst = 0 if sigmaSyst <= 0 else sigmaSyst ** 0.5
                hCFGenSyst.SetBinContent(iBin+1, muStat)
                hCFGenSyst.SetBinError(iBin+1, sigmaSyst)

            hCFGenTot.Write()
            hCFGenSyst.Write()

            gCFGenTot = ApplyCenterOfGravity(hCFGenTot, gGravities)
            gCFGenTot.SetName('gCFGenTot')
            gCFGenTot.Write()

            gCFGenSyst = ApplyCenterOfGravity(hCFGenSyst, gGravities)
            gCFGenSyst.SetName('gCFGenSyst')
            gCFGenSyst.Write()

            # Compute the average of the syst unc in k* in [0, 300] MeV
            systs = []
            for iBin in range(hCFGenSyst.FindBin(300)):
                systs.append(hCFGenSyst.GetBinError(iBin+1))
            systAvg = sum(systs) / len(systs)

            hCFGenSystAvg = hCFGenSyst.Clone('hCFGenSystAvg')
            for iBin in range(hCFGenSyst.FindBin(300)):
                hCFGenSystAvg.SetBinError(iBin + 1, systAvg)
            gCFGenSystAvg = ApplyCenterOfGravity(hCFGenSystAvg, gGravities)
            gCFGenSystAvg.SetName('gCFGenSystAvg')
            gCFGenSystAvg.Write()

            gBrackets = ComputeBinBrackets(hCFGenTot)
        gBrackets.Write()

    # Save purity
    # oFile.cd(comb)
    # for hPurity in dSEPurity:
    #     hPurity.Write()
    # for hPurity in dMEPurity:
    #     hPurity.Write()
        
    print(f'output saved in {oFileName}')
    oFile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarPi', 'DstarK'))
    parser.add_argument('--syst', action='store_true', default=False)
    parser.add_argument('--shift', action='store_true', default=False)
    parser.add_argument('--bs', type=int, default=0)
    args = parser.parse_args()

    ComputeGenCF(args)
