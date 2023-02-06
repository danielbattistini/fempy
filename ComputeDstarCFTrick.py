import argparse
import os
import ctypes
import numpy as np
import sys

from fempy.utils import Pair
from fempy import CorrelationFunction
import fempy

from ROOT import TFile, AliHFInvMassFitter, TDatabasePDG, TH1F, AliVertexingHFUtils, TH1D, TCanvas, gPad, gROOT, TF1, gInterpreter

gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter

def Gaus(x, pars):
    return pars[0] * np.e ** (-(x[0]-pars[1])**2/(2 * pars[2]**2))


def Pol0(x, pars):
    return pars[0]


def Pol1(x, pars):
    return pars[0] + pars[1]*x[0]


def GetNSigma(hist, func, name='nsigma'):
    hNSigma = hist.Clone(name)

    for iBin in range(hist.GetNbinsX()+1):
        iBinCenter = hist.GetXaxis().GetBinCenter(iBin+1)
        bc = hist.GetBinContent(iBin+1)
        bUnc = hist.GetBinError(iBin+1)
        expected = func.Eval(iBinCenter)
        hNSigma.SetBinContent(iBin+1, 0 if bUnc == 0 else (bc - expected) / bUnc)
    return hNSigma


def GetYieldsFromFitReweighted(hist, event):
    pass


# todo: to be removed in the future
# def GetYieldsFromFitUnreweighted(hist, event):
#     nKStarBins = round(hist.GetNbinsX() / kStarBW)
#     kStarEdges = [iKStar * kStarBW for iKStar in range(nKStarBins + 1)]
#     kStarMins = kStarEdges[:-1]
#     kStarMaxs = kStarEdges[1:]

#     nCanvases = (len(kStarMins) - 1) // nPadsPerCanvas + 1
#     cMasses = [TCanvas(f'c{event}Mass_{iC}', '', 1800, 1000) for iC in range(0, nCanvases)]
#     cResiduals = [TCanvas(f'c{event}Residuals_{iC}', '', 1800, 1000) for iC in range(0, nCanvases)]
#     cNSigmas = [TCanvas(f'c{event}NSigmas_{iC}', '', 1800, 1000) for iC in range(0, nCanvases)]

#     for canvases in [cMasses + cResiduals + cNSigmas]:
#         for canvas in canvases:
#             canvas = fempy.utils.DivideCanvas(canvas, nPadsPerCanvas)

#     currentDir = oFile.CurrentDirectory().load().GetName()

#     oFile.mkdir(f'{comb}/charm_mass')
#     oFile.cd(f'{comb}/charm_mass')

#     hYields = TH1D('hYields', ';#it{k}* (GeV/#it{c});Counts', nKStarBins, 0, kStarMaxs[-1])
#     hChi2 = TH1D('hChi2', ';#it{k}* (GeV/#it{c});#chi^{2}/NDF', nKStarBins, 0, kStarMaxs[-1])

#     fitters = []
#     for iKStarBin, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):

#         firstBin = hist.GetXaxis().FindBin(kStarMin * 1.0001)
#         lastBin = hist.GetXaxis().FindBin(kStarMax * 0.9999)

#         hCharmMass = hist.ProjectionY(f'h{event}CharmMass_kStar{kStarMin}_{kStarMax}', firstBin, lastBin)
#         hCharmMass.Write()

#         # bwOriginal = hist.GetYaxis().GetBinWidth(0) * 1000 #convert to MeV
#         charmMassRebin = round(charmMassBW / hist.GetYaxis().GetBinWidth(0) / 1000)

#         hMass = TH1F()

#         AliVertexingHFUtils.RebinHisto(hCharmMass, charmMassRebin).Copy(hMass)  # to cast TH1D to TH1F
#         hMass.SetTitle((f'{kStarMin:.0f} < #it{{k}}* < {kStarMax:.0f} MeV/#it{{c}};{massAxisTitle};'
#                         f'Counts per {charmMassBW:.1f} MeV/#it{{c}}^{{2}}'))

#         fitters.append(AliHFInvMassFitter(hMass, fitRange[0], fitRange[1], bkgFitFunc, AliHFInvMassFitter.kGaus))
#         fitters[-1].SetUseLikelihoodFit()
#         fitters[-1].SetInitialGaussianSigma(0.001)
#         fitters[-1].SetBoundGaussianMean(nominalMass, boundMassrange[0], boundMassrange[1])
#         status = fitters[-1].MassFitter(False)

#         if status != 1:  # fit failed
#             hYields.SetBinContent(iKStarBin+1, 0)
#             hYields.SetBinError(iKStarBin+1, 0)
#             hChi2.SetBinContent(iKStarBin+1, 0)
#             continue

#         # fMass = fitters[-1].GetMassFunc()
#         fSgn = fitters[-1].GetSignalFunc()
#         # fBkg = fitters[-1].GetBackgroundRecalcFunc()

#         hChi2.SetBinContent(iKStarBin + 1, fitters[-1].GetReducedChiSquare())

#         sgn, sgnUnc = ctypes.c_double(), ctypes.c_double()
#         fitters[-1].Signal(5, sgn, sgnUnc)

#         hYields.SetBinContent(iKStarBin + 1, 0 if fSgn == None else sgn)
#         hYields.SetBinError(iKStarBin + 1, 0 if fSgn == None else sgnUnc)

#         # Draw
#         cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
#         fitters[-1].DrawHere(gPad)

#         cResiduals[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
#         fitters[-1].DrawHistoMinusFit(gPad)

#         cNSigmas[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
#         hNSigmas = GetNSigma(hMass, fitters[-1].GetMassFunc(), f'h{event}NSigmas.pdf')
#         hNSigmas.GetYaxis().SetTitle('n_{#sigma}')
#         hNSigmas.GetYaxis().SetRangeUser(-5, 5)
#         hNSigmas.GetXaxis().SetRangeUser(fitRange[0], fitRange[1])
#         hNSigmas.DrawClone('hist')

#     for canvases in [cMasses, cResiduals, cNSigmas]:
#         for canvas in canvases:
#             canvas.Write()

#     for fitter in fitters:
#         fitter.Write()

#     oFile.cd(currentDir)
#     return hYields, hChi2


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', default="")
    parser.add_argument('inFileName')
    parser.add_argument('oDir')
    parser.add_argument('--suffix', default="")
    parser.add_argument('--kStarBW', default=50, type=float, help='units=MeV/c')
    parser.add_argument('--charmMassBW', default=0.5, type=float, help='units=MeV/c2')
    parser.add_argument('--fitMass', default=False, action='store_true')
    parser.add_argument('--aliFitter', default=False, action='store_true')
    parser.add_argument('--sgnFitFunc', default="")
    parser.add_argument('--bkgFitFunc', default="")

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
    nPadsPerCanvas = 6

    sgnFitFunc = args.sgnFitFunc
    if sgnFitFunc not in ["gaus", "hat"]:
        fempy.error("signal fit function not implemented")

    bkgFitFunc = args.bkgFitFunc
    if bkgFitFunc not in ["powex"]:
        fempy.error("background fit function not implemented")

    if pair.name == 'DstarPi':
        nominalMass = TDatabasePDG.Instance().GetParticle(413).Mass() - TDatabasePDG.Instance().GetParticle(421).Mass()
        massAxisTitle = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
        boundMassrange = [0.144, 0.146]
        fitRange = [0.141, 0.154]
    elif pair.name == 'DPi':
        nominalMass = TDatabasePDG.Instance().GetParticle(411).Mass()
        massAxisTitle = pair.heavy_mass_label
        boundMassrange = [1.85, 1.9]
        fitRange = [1.8, 1.95]
    else:
        fempy.error('hadron not implemented')

    for comb in combs:
        oFile.mkdir(f'{comb}/fits')

        if args.fitMass:  # todo: adapt for Dmeson-pi Dmeson-K
            hDistrPurityFromSgnFuncInt = {}
            hDistrPurityFromDataMinusBkg = {}
            for event in ['SE', 'ME']:
                oFile.cd(f'{comb}')

                hMasses = [h for h in fempy.utils.GetObjsInDir(inFile.Get(f'{comb}/{event}')) if 'hCharmMass0' in h.GetName()]

                nKStarBins = len(hMasses)
                kStarBW = round(3000/nKStarBins)

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

                hPurityFromSgnFuncInt = TH1D(f'h{event}PurityFromSgnFuncInt', ';#it{k}* (GeV/#it{c});Purity', nKStarBins, 0, kStarMaxs[-1])
                hPurityFromDataMinusBkg = TH1D(f'h{event}PurityFromDataMinusBkg', ';#it{k}* (GeV/#it{c});Purity', nKStarBins, 0, kStarMaxs[-1])
                hBkg = TH1D(f'h{event}Bkg', ';#it{k}* (GeV/#it{c});Counts', nKStarBins, 0, kStarMaxs[-1])
                hYieldsFromSgnFuncInt = TH1D(f'h{event}YieldsFromFit', ';#it{k}* (GeV/#it{c});Counts', nKStarBins, 0, kStarMaxs[-1])
                hYieldsFromDataMinusBkg = TH1D(f'h{event}YieldsFromDataMinusBkg', ';#it{k}* (GeV/#it{c});Counts', nKStarBins, 0, kStarMaxs[-1])
                hWidth = TH1D(f'h{event}Width', ';#it{k}* (GeV/#it{c});Width (GeV/#it{c})', nKStarBins, 0, kStarMaxs[-1])
                hChi2 = TH1D(f'h{event}Chi2', ';#it{k}* (GeV/#it{c});#chi^{2}/NDF', nKStarBins, 0, kStarMaxs[-1])

                fitters = []

                for iKStarBin, (hCharmMass, kStarMin, kStarMax) in enumerate(zip(hMasses, kStarMins, kStarMaxs)):
                    print(f"\n\nProcessing {comb} {event} kStar: [{kStarMin:.0f}, {kStarMax:.0f}]")

                    charmMassRebin = round(charmMassBW / hCharmMass.GetXaxis().GetBinWidth(0) / 1000)
                    nSigma = 2

                    if args.aliFitter:
                        hMass = TH1F()

                        AliVertexingHFUtils.RebinHisto(hCharmMass, charmMassRebin).Copy(hMass)  # to cast TH1D to TH1F
                        hMass.SetTitle((f'{kStarMin:.0f} < #it{{k}}* < {kStarMax:.0f} MeV/#it{{c}};{massAxisTitle};'
                                        f'Counts per {charmMassBW:.1f} MeV/#it{{c}}^{{2}}'))

                        if bkgFitFunc == "powex":
                            aliBkgFitFunc = AliHFInvMassFitter.kPowEx
                        if sgnFitFunc == "gaus":
                            aliSgnFitFunc = AliHFInvMassFitter.kGaus
                        elif sgnFitFunc == "hat":
                            aliSgnFitFunc = AliHFInvMassFitter.k2GausSigmaRatioPar

                        fitters.append(AliHFInvMassFitter(hMass, fitRange[0]*1.0001, fitRange[1]*0.9999, aliBkgFitFunc, aliSgnFitFunc))
                        fitters[-1].SetUseLikelihoodFit()
                        fitters[-1].SetInitialGaussianSigma(0.0006)
                        fitters[-1].SetBoundGaussianMean(nominalMass, boundMassrange[0], boundMassrange[1])
                        # if aliSgnFitFunc == AliHFInvMassFitter.k2GausSigmaRatioPar:
                        #     # fitters[-1].SetFixGaussianMean(0.145)
                        #     # fitters[-1].SetFixGaussianSigma(0.0006)
                        #     # fitters[-1].SetInitialRatio2GausSigma(0.7)

                        #     fitters[-1].SetFixRatio2GausSigma(0.6)
                        # fitters[-1].SetFixRatio2GausSigma(0.1)
                        # fitters[-1].SetBoundGaussianSigma(0.0005, 0.001)
                        status = fitters[-1].MassFitter(False)
                        # ccc = TCanvas("lallero", "Llall", 600, 600)
                        # ccc.SaveAs("test.png")

                        if status != 1:  # fit failed
                            hYieldsFromSgnFuncInt.SetBinContent(iKStarBin+1, 0)
                            hYieldsFromSgnFuncInt.SetBinError(iKStarBin+1, 0)
                            hPurityFromSgnFuncInt.SetBinError(iKStarBin+1, 0)
                            hPurityFromSgnFuncInt.SetBinError(iKStarBin+1, 0)
                            hPurityFromDataMinusBkg.SetBinError(iKStarBin+1, 0)
                            hPurityFromDataMinusBkg.SetBinError(iKStarBin+1, 0)
                            hChi2.SetBinContent(iKStarBin+1, 0)
                            cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
                            hCharmMass.DrawClone()
                            continue
                        chi2 = fitters[-1].GetReducedChiSquare()
                        fMass = fitters[-1].GetMassFunc()
                        fSgn = fitters[-1].GetSignalFunc()
                        fBkg = fitters[-1].GetBackgroundRecalcFunc()
                        sgn, sgnUnc = ctypes.c_double(), ctypes.c_double()
                        bkg, bkgUnc = ctypes.c_double(), ctypes.c_double()
                        fitters[-1].Signal(2, sgn, sgnUnc)
                        width = fitters[-1].GetSigma()
                        widthUnc = fitters[-1].GetSigmaUncertainty()
                        fitters[-1].Background(2, bkg, bkgUnc)

                        cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
                        fitters[-1].DrawHere(gPad, 2)

                        cResiduals[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
                        fitters[-1].DrawHistoMinusFit(gPad)

                        sgnFromSgnFuncInt = sgn.value
                        sgnFromSgnFuncIntUnc = sgnUnc.value
                        bkg = bkg.value
                        bkgUnc = bkgUnc.value
                        

                    else:
                        hCharmMass.SetTitle(f'{kStarMin:.0f}  <  #it{{k}}* < {kStarMax:.0f} MeV/#it{{c}}')
                        charmMassRebin = round(charmMassBW / hCharmMass.GetXaxis().GetBinWidth(0) / 1000)
                        hCharmMass.Rebin(charmMassRebin)
                        hCharmMass.GetXaxis().SetRangeUser(0.135, 0.160)
                        fitters.append(MassFitter(hCharmMass, sgnFitFunc, "powex", fitRange[0], fitRange[1]))
                        status = fitters[-1].Fit()


                        chi2 = fitters[-1].GetChi2Ndf()
                        sgnFromSgnFuncInt = fitters[-1].GetSignal(nSigma, 'sgn_int')
                        sgnFromSgnFuncIntUnc = fitters[-1].GetSignalUnc(nSigma, 'sgn_int')
                        sgnFromDataMinusBkgFuncInt = fitters[-1].GetSignal(nSigma, 'data_minus_bkg')
                        sgnFromDataMinusBkgFuncIntUnc = fitters[-1].GetSignalUnc(nSigma, 'data_minus_bkg')
                        
                        bkg = fitters[-1].GetBackground(nSigma)
                        bkgUnc = fitters[-1].GetBackgroundUnc(nSigma)

                        width = fitters[-1].GetWidth()
                        widthUnc = fitters[-1].GetWidthUnc()

                        # Draw
                        cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
                        fitters[-1].Draw(gPad, 'data_minus_bkg')
                        
                    hChi2.SetBinContent(iKStarBin + 1, chi2)
                    hBkg.SetBinContent(iKStarBin + 1, bkg)
                    hBkg.SetBinError(iKStarBin + 1, bkgUnc)
                    hYieldsFromSgnFuncInt.SetBinContent(iKStarBin + 1, sgnFromSgnFuncInt)
                    hYieldsFromSgnFuncInt.SetBinError(iKStarBin + 1, sgnFromSgnFuncIntUnc)

                    firstBin = hCharmMass.GetXaxis().FindBin((nominalMass - nSigma * width)*1.0001)
                    lastBin = hCharmMass.GetXaxis().FindBin((nominalMass + nSigma * width)*0.9999)
                    sgnFromHist = hCharmMass.Integral(firstBin, lastBin)
                    sgnFromHistUnc = np.sqrt(sum((hCharmMass.GetBinContent(bin) for bin in range(firstBin, lastBin+1))))

                    hYieldsFromDataMinusBkg.SetBinContent(iKStarBin + 1, sgnFromHist - bkg)
                    hYieldsFromDataMinusBkg.SetBinError(iKStarBin + 1, np.sqrt(sgnFromHistUnc**2 + bkgUnc**2))

                    hWidth.SetBinContent(iKStarBin + 1, width)
                    hWidth.SetBinError(iKStarBin + 1, widthUnc)
                    hPurityFromSgnFuncInt.SetBinContent(iKStarBin + 1, sgnFromSgnFuncInt/(sgnFromSgnFuncInt + bkg))
                    hPurityFromSgnFuncInt.SetBinError(iKStarBin + 1, np.sqrt(bkg**2 * sgnFromSgnFuncIntUnc**2 + sgnFromSgnFuncInt**2 * bkgUnc ** 2) / (sgnFromSgnFuncInt + bkg)**2)
                    if not args.aliFitter:
                        print("akjnkajnk", sgnFromDataMinusBkgFuncInt)
                        hPurityFromDataMinusBkg.SetBinContent(iKStarBin + 1, sgnFromDataMinusBkgFuncInt/(sgnFromDataMinusBkgFuncInt + bkg))
                        hPurityFromDataMinusBkg.SetBinError(iKStarBin + 1, np.sqrt(bkg**2 * sgnFromDataMinusBkgFuncIntUnc**2 + sgnFromDataMinusBkgFuncInt**2 * bkgUnc ** 2) / (sgnFromDataMinusBkgFuncInt + bkg)**2)

                hBkg.Write()
                hYieldsFromSgnFuncInt.Write()
                hPurityFromSgnFuncInt.Write()
                if not args.aliFitter:
                    hYieldsFromDataMinusBkg.Write()
                    hPurityFromDataMinusBkg.Write()

                hWidth.Write()
                hChi2.Write()

                hDistrPurityFromSgnFuncInt[event] = inFile.Get(f'{comb}/{event}/sgn/hCharmMassVsKStar0').ProjectionX()
                hDistrPurityFromSgnFuncInt[event].SetName(f'h{event}PurityRewFromSgnFuncInt')
                hDistrPurityFromSgnFuncInt[event].Rebin(round(kStarBW/hDistrPurityFromSgnFuncInt[event].GetXaxis().GetBinWidth(1)))
                for iKStarBin in range(nKStarBins):
                    hDistrPurityFromSgnFuncInt[event].SetBinContent(iKStarBin, hDistrPurityFromSgnFuncInt[event].GetBinContent(iKStarBin + 1) * hPurityFromSgnFuncInt.GetBinContent(iKStarBin+1))
                hDistrPurityFromSgnFuncInt[event].Write()

                if not args.aliFitter:
                    hDistrPurityFromDataMinusBkg[event] = inFile.Get(f'{comb}/{event}/sgn/hCharmMassVsKStar0').ProjectionX()
                    hDistrPurityFromDataMinusBkg[event].SetName(f'h{event}PurityRewFromDataMinusBkg')
                    hDistrPurityFromDataMinusBkg[event].Rebin(round(kStarBW/hDistrPurityFromDataMinusBkg[event].GetXaxis().GetBinWidth(1)))
                    for iKStarBin in range(nKStarBins):
                        hDistrPurityFromDataMinusBkg[event].SetBinContent(iKStarBin, hDistrPurityFromDataMinusBkg[event].GetBinContent(iKStarBin + 1) * hPurityFromSgnFuncInt.GetBinContent(iKStarBin+1))
                    hDistrPurityFromDataMinusBkg[event].Write()

                oFile.cd(f'{comb}/fits')
                oInvMassFitsFileName = os.path.join(args.oDir, f"InvMassFits_{comb}_{event}.pdf" if args.suffix == '' else f'InvMassFits_{comb}_{event}_{args.suffix}.pdf')
                cMasses[0].SaveAs(f'{oInvMassFitsFileName}[')
                for canvas in cMasses:
                    canvas.Write()
                    canvas.SaveAs(f'{oInvMassFitsFileName}')
                cMasses[-1].SaveAs(f'{oInvMassFitsFileName}]')

                if args.aliFitter:
                    for canvas in cResiduals:
                        canvas.Write()
                oFile.cd(f'{comb}')

            hCF = CorrelationFunction(se=hDistrPurityFromSgnFuncInt['SE'], me=hDistrPurityFromSgnFuncInt['ME'], norm=pair.norm_range).get_cf()
            hCF.Write('hCFPurityRewFromSgnFuncInt')
            if not args.aliFitter:
                hCF = CorrelationFunction(se=hDistrPurityFromDataMinusBkg['SE'], me=hDistrPurityFromDataMinusBkg['ME'], norm=pair.norm_range).get_cf()
                hCF.Write('hCFPurityRewFromDataMinusBkg')

        for region in regions:
            hSEMultVsKStar = inFile.Get(f'{comb}/SE/{region}/hCharmMassVsKStar0')
            hSEMultVsKStar.SetName('hSECharmMassVsKStar0')
            hMEMultVsKStar = inFile.Get(f'{comb}/ME/{region}/hCharmMassVsKStar0')
            hMEMultVsKStar.SetName('hMECharmMassVsKStar0')

            hSE = hSEMultVsKStar.ProjectionX()
            hSE.Write(f'hSEstd_{region}')
            hME = hMEMultVsKStar.ProjectionX()
            hME.Write(f'hMEstd_{region}')

            hSE.Rebin(round(kStarBW / hSE.GetXaxis().GetBinWidth(0)))
            hME.Rebin(round(kStarBW / hME.GetXaxis().GetBinWidth(0)))

            hCF = CorrelationFunction(se=hSE, me=hME, norm=pair.norm_range, units='MeV').get_cf()
            hCF.Write(f'hCFstd_{region}')
    oFile.Close()
    print(f"output saved in {oFileName}")
