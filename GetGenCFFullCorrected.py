import sys
import yaml

import argparse

from ROOT import TFile, TF1, TGraphErrors, TDatabasePDG, gInterpreter, TCanvas, kBlue, TLatex, gStyle, TLegend
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


def LoadGravities(hME, kStarBW):
    nKStarBins = round(hME.GetNbinsX()/kStarBW)
    xCF = [Average(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
    xCFUnc = [StdDev(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
    yCF = [1 for _ in range(nKStarBins)]
    yCFUnc = [0 for _ in range(nKStarBins)]
    gGravities = TGraphErrors(1)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarPi', 'DstarK'))
    parser.add_argument('--syst', action='store_true', default=False)
    parser.add_argument('--stat', action='store_true', default=False)
    args = parser.parse_args()

    gStyle.SetOptStat(0)

    kStarBW = 50  # MeV/c
    if args.pair == 'DstarK':
        fitRanges = [[10, 450], [10, 400], [10, 500]]
        inFileData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq.root'
        config = '/home/daniel/an/DstarK/cfg_gencf_DstarK_50MeV.yml'

        lightMass = TDatabasePDG.Instance().GetParticle(321).Mass()

        weights1 = [0.78, 0.80, 0.77]
        radii1 = [0.86, 0.95, 0.79]
        radii2 = [2.03, 2.22, 1.91]
    else:
        print("not implemented")
        sys.exit()

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
        hSEMC.Write()
        hMEMC.Write()
        hCFMC.SetName('hCFMC')
        hCFMC.Write()

        regions = fempy.utils.GetRegions(inFileData.Get('sc/SE'))
        dCFData = {}
        for region in regions:
            hSEData = inFileData.Get(f'{comb}/SE/{region}/hCharmMassVsKStar0').ProjectionX(f'hSE_{region}')
            hMEData = inFileData.Get(f'{comb}/ME/{region}/hCharmMassVsKStar0').ProjectionX(f'hME_{region}')

            rebinFactor = round(kStarBW/hSEData.GetBinWidth(1))
            hSEData.Rebin(rebinFactor)
            hMEData.Rebin(rebinFactor)
            hSEData.Scale(ComputeNormFactor(hSEData, hMEData, normRange[0], normRange[1]))
            hCFData = hSEData/hMEData
            hSEData.Write()
            hMEData.Write()
            hCFData.SetName(f'hCF_{region}')
            hCFData.Write()
            dCFData[region] = hCFData

        # Prefit the background
        hCFNorm = (dCFData['sgn'] - lamPar['sb'] * dCFData['sbr']) / (hCFMC * (lamPar['gen'] + lamPar['flat']))
        hCFNorm.SetName('hCFNorm')
        hCFNorm.Write()

        # Compute center of gravity of the bins in the ME
        hGravities = inFileData.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hGravities')
        print(kStarBW)
        gGravities = LoadGravities(hGravities, kStarBW)

        # Fit the baseline
        fBaseLine = TF1('fBaseLine', '[0]', 0, 3000)
        ApplyCenterOfGravity(hCFNorm, gGravities).Fit(fBaseLine, 'Q', '', 300, 1000)
        blNorm = fBaseLine.GetParameter(0)

        # Compute the total background mode
        hCFBkg = lamPar['sb'] * dCFData['sbr'] + blNorm * hCFMC * (lamPar['gen'] + lamPar['flat'])
        hCFBkg.SetName('hCFBkg')
        hCFBkg.Write()

        # Compute the Gen CF
        hCFFlat = MakeFlatHist(dCFData['sgn'], lamPar['flat'], 'hCFFlat')
        hCFGen = (dCFData['sgn'] - lamPar['sb'] * dCFData['sbr'])/(blNorm * hCFMC) - hCFFlat
        hCFGen.Scale(1./lamPar['gen'])
        hCFGen.SetName('hCFGen')
        hCFGen.Write()

        # Compute the coulomb-only CF with Lednicky
        fitRange = fitRanges[0]
        fCoulomb = TF1("fCoulomb", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
        fCoulomb.FixParameter(0, radii1[0])
        fCoulomb.FixParameter(1, radii2[0])
        fCoulomb.FixParameter(2, weights1[0])
        fCoulomb.FixParameter(3, 0)
        fCoulomb.FixParameter(4, 0.)
        fCoulomb.FixParameter(5, 0.)
        fCoulomb.FixParameter(6, 1 if comb == 'sc' else -1)
        fCoulomb.FixParameter(7, RedMass(heavyMass, lightMass) * 1000)
        fCoulomb.SetLineColor(kBlue)

        # Fit the CF
        gCFGen = ApplyCenterOfGravity(hCFGen, gGravities)
        fWeightedLL = TF1("fWeightedLL", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
        fWeightedLL.FixParameter(0, radii1[0])
        fWeightedLL.FixParameter(1, radii2[0])
        fWeightedLL.FixParameter(2, weights1[0])
        fWeightedLL.SetParameter(3, 0.1)
        fWeightedLL.SetParLimits(3, -1, 1)
        fWeightedLL.FixParameter(4, 0.)
        fWeightedLL.FixParameter(5, 0.)
        fWeightedLL.FixParameter(6, 1 if comb == 'sc' else -1)
        fWeightedLL.FixParameter(7, RedMass(heavyMass, lightMass) * 1000)

        status = gCFGen.Fit(fWeightedLL, "SMRQ+0").Status()
        chi2ndf = fWeightedLL.GetChisquare()/fWeightedLL.GetNDF()
        scattLen = fWeightedLL.GetParameter(3)
        scattLenUnc = fWeightedLL.GetParError(3)

        cFit = TCanvas(f'cFit_{comb}', '', 600, 600)
        cFit.DrawFrame(0, 0, 500, 2)
        hCFGen.Draw('pe')
        fWeightedLL.Draw('same')
        fCoulomb.Draw('same')

        step = 0.05
        iVar = 0
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
        cFit.Write()

    print(f'output saved in {oFileName}')
    oFile.Close()