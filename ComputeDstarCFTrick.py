import argparse
import os
import ctypes

from ROOT import TFile, AliHFInvMassFitter, TDatabasePDG, TH1F, AliVertexingHFUtils, TH1D, TCanvas, gPad

import fempy
from fempy import CorrelationFunction
from fempy.utils import Pair


def GetNSigma(hist, func, name='nsigma'):
    hNSigma = hist.Clone(name)

    for iBin in range(hist.GetNbinsX()+1):
        iBinCenter = hist.GetXaxis().GetBinCenter(iBin+1)
        bc = hist.GetBinContent(iBin+1)
        bUnc = hist.GetBinError(iBin+1)
        expected = func.Eval(iBinCenter)
        hNSigma.SetBinContent(iBin+1, 0 if bUnc == 0 else (bc - expected) / bUnc)
    return hNSigma


def GetYieldsFromFit(hist, event):
    nKStarBins = round(hist.GetNbinsX() / kStarBW)
    kStarEdges = [iKStar * kStarBW for iKStar in range(nKStarBins + 1)]
    kStarMins = kStarEdges[:-1]
    kStarMaxs = kStarEdges[1:]

    nCanvases = (len(kStarMins) - 1) // nPadsPerCanvas + 1
    cMasses = [TCanvas(f'c{event}Mass_{iC}', '', 1920, 1200) for iC in range(0, nCanvases)]
    cResiduals = [TCanvas(f'c{event}Residuals_{iC}', '', 1920, 1200) for iC in range(0, nCanvases)]
    cNSigmas = [TCanvas(f'c{event}NSigmas_{iC}', '', 1920, 1200) for iC in range(0, nCanvases)]

    for canvases in [cMasses + cResiduals + cNSigmas]:
        for canvas in canvases:
            canvas = fempy.utils.DivideCanvas(canvas, nPadsPerCanvas)

    currentDir = oFile.CurrentDirectory().load().GetName()

    oFile.mkdir(f'{comb}/charm_mass')
    oFile.cd(f'{comb}/charm_mass')

    hYields = TH1D('hYields', ';#it{k} (GeV/#it{c}^{2});Counts', nKStarBins, 0, kStarMaxs[-1])
    hChi2 = TH1D('hChi2', ';#it{k} (GeV/#it{c}^{2});#chi^{2}/NDF', nKStarBins, 0, kStarMaxs[-1])

    fitters = []
    for iKStarBin, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):

        firstBin = hist.GetXaxis().FindBin(kStarMin * 1.0001)
        lastBin = hist.GetXaxis().FindBin(kStarMax * 0.9999)

        hCharmMass = hist.ProjectionY(f'h{event}CharmMass_kStar{kStarMin}_{kStarMax}', firstBin, lastBin)
        hCharmMass.Write()

        # bwOriginal = hist.GetYaxis().GetBinWidth(0) * 1000 #convert to MeV
        charmMassRebin = round(charmMassBW / hist.GetYaxis().GetBinWidth(0) / 1000)

        hMass = TH1F()

        AliVertexingHFUtils.RebinHisto(hCharmMass, charmMassRebin).Copy(hMass)  # to cast TH1D to TH1F
        hMass.SetTitle((f'{kStarMin:.0f} < #it{{k}}* < {kStarMax:.0f} MeV/#it{{c}};{massAxisTitle};'
                        f'Counts per {charmMassBW:.1f} MeV/#it{{c}}^{{2}}'))

        fitters.append(AliHFInvMassFitter(hMass, fitRange[0], fitRange[1], bkgFitFunc, AliHFInvMassFitter.kGaus))
        fitters[-1].SetUseLikelihoodFit()
        fitters[-1].SetInitialGaussianSigma(0.001)
        fitters[-1].SetBoundGaussianMean(nominalMass, boundMassrange[0], boundMassrange[1])
        status = fitters[-1].MassFitter(False)

        if status != 1:  # fit failed
            hYields.SetBinContent(iKStarBin+1, 0)
            hYields.SetBinError(iKStarBin+1, 0)
            hChi2.SetBinContent(iKStarBin+1, 0)
            continue

        # fMass = fitters[-1].GetMassFunc()
        fSgn = fitters[-1].GetSignalFunc()
        # fBkg = fitters[-1].GetBackgroundRecalcFunc()

        hChi2.SetBinContent(iKStarBin + 1, fitters[-1].GetReducedChiSquare())

        sgn, sgnUnc = ctypes.c_double(), ctypes.c_double()
        fitters[-1].Signal(5, sgn, sgnUnc)

        hYields.SetBinContent(iKStarBin + 1, 0 if fSgn == None else sgn)
        hYields.SetBinError(iKStarBin + 1, 0 if fSgn == None else sgnUnc)

        # Draw
        cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
        fitters[-1].DrawHere(gPad)

        cResiduals[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
        fitters[-1].DrawHistoMinusFit(gPad)

        cNSigmas[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
        hNSigmas = GetNSigma(hMass, fitters[-1].GetMassFunc(), f'h{event}NSigmas.pdf')
        hNSigmas.GetYaxis().SetTitle('n_{#sigma}')
        hNSigmas.GetYaxis().SetRangeUser(-5, 5)
        hNSigmas.GetXaxis().SetRangeUser(fitRange[0], fitRange[1])
        hNSigmas.DrawClone('hist')

    for canvases in [cMasses, cResiduals, cNSigmas]:
        for canvas in canvases:
            canvas.Write()

    for fitter in fitters:
        fitter.Write()

    oFile.cd(currentDir)
    return hYields, hChi2


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', default="")
    parser.add_argument('inFileName')
    parser.add_argument('oDir')
    parser.add_argument('--suffix', default="")
    parser.add_argument('--kStarBW', default=50, type=float, help='units=MeV/c')
    parser.add_argument('--charmMassBW', default=1, type=float, help='units=MeV/c2')
    parser.add_argument('--fitMass', default=False, action='store_true')

    args = parser.parse_args()

    inFile = TFile(args.inFileName)
    oFileName = os.path.join(args.oDir, "RawCF.root" if args.suffix == '' else f'RawCF_{args.suffix}.root')
    oFile = TFile(oFileName, 'recreate')

    combs = fempy.utils.io.GetSubdirsInDir(inFile)
    combs = [comb for comb in combs if len(comb) == 2]
    regions = fempy.utils.io.GetSubdirsInDir(inFile.Get(f'{combs[0]}/SE'))
    pair = Pair(args.pair)
    kStarBW = args.kStarBW
    charmMassBW = args.charmMassBW
    nPadsPerCanvas = 12

    if pair.name == 'DstarPi':
        nominalMass = TDatabasePDG.Instance().GetParticle(413).Mass() - TDatabasePDG.Instance().GetParticle(421).Mass()
        massAxisTitle = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
        boundMassrange = [0.144, 0.146]
        fitRange = [0.13, 0.16]
        bkgFitFunc = AliHFInvMassFitter.kPowEx
    elif pair.name == 'DPi':
        nominalMass = TDatabasePDG.Instance().GetParticle(411).Mass()
        print(nominalMass)
        massAxisTitle = pair.heavy_mass_label
        boundMassrange = [1.85, 1.9]
        fitRange = [1.8, 1.95]
        bkgFitFunc = AliHFInvMassFitter.kPol2
    else:
        fempy.error('hadron not implemented')

    for comb in combs:
        oFile.mkdir(comb)
        oFile.cd(comb)

        if args.fitMass:  # todo: adapt for Dmeson-pi Dmeson-K
            hSEMultVsKStar = inFile.Get(f'{comb}/SE/hCharmMassVsKStar0')
            inFile.ls()
            hSEMultVsKStar.SetName('hSEMultVsKStar')
            hSEYieldSignal, hSEChi2 = GetYieldsFromFit(hSEMultVsKStar, 'SE')

            hMEMultVsKStar = inFile.Get(f'{comb}/ME/hCharmMassVsKStar0')
            hMEMultVsKStar.SetName('hMEMultVsKStar')
            hMEYieldSignal, hMEChi2 = GetYieldsFromFit(hMEMultVsKStar, 'ME')

            hSEYieldSignal.Write('hSEYields')
            hSEChi2.Write('hSEChi2')

            hMEYieldSignal.Write('hMEYields')
            hMEChi2.Write('hMEChi2')

            hCF = CorrelationFunction(se=hSEYieldSignal, me=hMEYieldSignal, norm=pair.norm_range).get_cf()
            hCF.Write('hCFfit')

        for region in regions:
            hSEMultVsKStar = inFile.Get(f'{comb}/SE/{region}/hCharmMassVsKStar0')
            hSEMultVsKStar.SetName('hSECharmMassVsKStar0')
            hMEMultVsKStar = inFile.Get(f'{comb}/ME/{region}/hCharmMassVsKStar0')
            hMEMultVsKStar.SetName('hMECharmMassVsKStar0')

            hSE = hSEMultVsKStar.ProjectionX()
            hSE.Write('SE')
            hME = hMEMultVsKStar.ProjectionX()
            hME.Write('ME')

            hSE.Rebin(round(kStarBW / hSE.GetXaxis().GetBinWidth(0)))
            hME.Rebin(round(kStarBW / hME.GetXaxis().GetBinWidth(0)))

            hCF = CorrelationFunction(se=hSE, me=hME, norm=pair.norm_range, units='MeV').get_cf()
            hCF.Write(f'hCFstd_{region}')

    oFile.Close()
    print(f"output saved in {oFileName}")
