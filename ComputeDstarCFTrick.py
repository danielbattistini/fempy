import argparse
import os
import ctypes
import numpy as np
import sys


from ROOT import TFile, AliHFInvMassFitter, TDatabasePDG, TH1F, AliVertexingHFUtils, TH1D, TCanvas, gPad, gROOT, TF1, gInterpreter

gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')

from ROOT import MassFitter

# a = MassFitter()
# sys.exit()


import fempy
from fempy import CorrelationFunction
from fempy.utils import Pair


def Gaus(x, pars):
    return pars[0] * np.e ** (-(x[0]-pars[1])**2/(2 * pars[2]**2))

def Pol0(x, pars):
    return pars[0]

def Pol1(x, pars):
    print(type(pars))
    return pars[0] + pars[1]*x[0]
# class MassFitter:
#     def __init__(self, hist, sgnFuncName, bkgFuncName, fitRange):
#         sgnFuncs = {
#             'gaus' : Gaus
#         }
#         self.hist = hist
#         self.sgnFunc = sgnFuncs[sgnFuncName]
#         self.bkgFunc = bkgFuncName
#         self.fitRange = fitRange
#         self.fSgn = None
#         self.fBkg = None
#         self.fTot = None
#         self.TotFit = None


    # def MakeTotFitFunc(self):
    #     sgnFunc.Eval()




    # def Fit(self):
    #     print(self.sgnFunc)
    #     self.fSgn = TF1("fSgn", Gaus, fitRange[0], fitRange[1], 3)
        
    #     self.fBkg = TF1("fBkg", Pol1, fitRange[0], fitRange[1], 2)
        

    #     nSigPars = self.fSgn.GetNpar()
    #     nBkgPars = self.fBkg.GetNpar()
    #     nPars = nSigPars + nBkgPars

    #     # self.TotFit = lambda x, pars : self.fSgn.EvalPar(x, pars) + self.fBkg.EvalPar(x, np.asarray(pars[:], 'd'))
    #     self.TotFit = lambda x, pars : Gaus(x, pars) + Pol1(x, pars[:])

    #     self.fTot = TF1("fTot", self.TotFit, fitRange[0], fitRange[1], nPars)
    #     self.fTot.SetParameter(0, 60)
    #     # self.fTot.SetParLimits(0, 0, 1000)
    #     self.fTot.SetParameter(1, 0.145)
    #     self.fTot.SetParLimits(1, 0.144, 0.146)
    #     self.fTot.SetParameter(2, 0.001)
    #     self.fTot.SetParLimits(2, 0.0002, 0.002)
    #     self.fTot.SetParameter(3, 20)
    #     self.fTot.SetParLimits(3, 0, 30)
    #     self.fTot.SetParameter(4, 1000)
    #     self.fTot.SetParLimits(4, 0, 3000)
        
    #     self.hist.Fit(self.fTot, "MR+", "", fitRange[0], fitRange[1])



    #     # fBkg = TF1('fBkg', "pol0", fitRange[0], fitRange[1], 1)
    #     # fBkg.SetParLimits(0, 0, 1000000)

    #     # def TotFitt(self, x):
    #     #     return self.fSgn.Eval(x[0])
    #     # nPars = 3
    #     # # def TotFit(x,  pars):
    #     # #     return self.fSgn.Eval(x[0])
        
    #     # self.fTot = TF1('fTot', TotFitt, fitRange[0], fitRange[1], nPars)

    #     # self.hist.Fit(self.fTot, "", "", fitRange[0], fitRange[1])

    # def Draw(self):
    #     self.hist.Draw()


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
def GetYieldsFromFitUnreweighted(hist, event):

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
    parser.add_argument('--charmMassBW', default=0.5, type=float, help='units=MeV/c2')
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
        fitRange = [0.140, 0.16]
        sgnFitFunc = AliHFInvMassFitter.kGaus
        # sgnFitFunc = AliHFInvMassFitter.k2GausSigmaRatioPar
        bkgFitFunc = AliHFInvMassFitter.kPowEx
    elif pair.name == 'DPi':
        nominalMass = TDatabasePDG.Instance().GetParticle(411).Mass()
        print(nominalMass)
        massAxisTitle = pair.heavy_mass_label
        boundMassrange = [1.85, 1.9]
        fitRange = [1.8, 1.95]
        sgnFitFunc = AliHFInvMassFitter.kGaus
        bkgFitFunc = AliHFInvMassFitter.kPol2
    else:
        fempy.error('hadron not implemented')

    for comb in combs:
        oFile.mkdir(comb)
        oFile.cd(comb)

        if args.fitMass:  # todo: adapt for Dmeson-pi Dmeson-K
            hDistr = {}
            for event in ['SE', 'ME']:
                hMasses = [h for h in fempy.utils.GetObjsInDir(inFile.Get(f'{comb}/{event}')) if 'hCharmMass0' in h.GetName()]
                print(len(hMasses))
                
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

                oFile.mkdir(f'{comb}')
                oFile.cd(f'{comb}')

                hPurity = TH1D(f'h{event}Purity', ';#it{k} (GeV/#it{c}^{2});Counts', nKStarBins, 0, kStarMaxs[-1])
                hYields = TH1D(f'h{event}Yields', ';#it{k} (GeV/#it{c}^{2});Counts', nKStarBins, 0, kStarMaxs[-1])
                hChi2 = TH1D(f'h{event}Chi2', ';#it{k} (GeV/#it{c}^{2});#chi^{2}/NDF', nKStarBins, 0, kStarMaxs[-1])

                fitters = []
                # gROOT.SetBatch(False)
                for iKStarBin, (hCharmMass, kStarMin, kStarMax) in enumerate(zip(hMasses, kStarMins, kStarMaxs)):
                    # charmMassRebin = round(charmMassBW / hCharmMass.GetXaxis().GetBinWidth(0) / 1000)
                    # hCharmMass.Rebin(charmMassRebin)
                    # hCharmMass.GetXaxis().SetRangeUser(0.135, 0.160)
                    # fitter = MassFitter(hCharmMass, "hat", "powex", fitRange[0], fitRange[1])
                    # fitter.Fit()
                    # ccc = TCanvas("lallero", "Llall", 600, 600)
                    # fitter.Draw()
                    # ccc.SaveAs("test.png")
                    break

                for iKStarBin, (hCharmMass, kStarMin, kStarMax) in enumerate(zip(hMasses, kStarMins, kStarMaxs)):
                    charmMassRebin = round(charmMassBW / hCharmMass.GetXaxis().GetBinWidth(0) / 1000)
                    print(charmMassBW, hCharmMass.GetXaxis().GetBinWidth(0), charmMassRebin)
                    hMass = TH1F()

                    AliVertexingHFUtils.RebinHisto(hCharmMass, charmMassRebin).Copy(hMass)  # to cast TH1D to TH1F
                    hMass.SetTitle((f'{kStarMin:.0f} < #it{{k}}* < {kStarMax:.0f} MeV/#it{{c}};{massAxisTitle};'
                                    f'Counts per {charmMassBW:.1f} MeV/#it{{c}}^{{2}}'))

                    fitters.append(AliHFInvMassFitter(hMass, fitRange[0], fitRange[1], bkgFitFunc, sgnFitFunc))
                    fitters[-1].SetUseLikelihoodFit()
                    fitters[-1].SetInitialGaussianSigma(0.0006)
                    fitters[-1].SetBoundGaussianMean(nominalMass, boundMassrange[0], boundMassrange[1])
                    if sgnFitFunc == AliHFInvMassFitter.k2GausSigmaRatioPar:
                        # fitters[-1].SetFixGaussianMean(0.145)
                        # fitters[-1].SetFixGaussianSigma(0.0006)
                        # fitters[-1].SetInitialRatio2GausSigma(0.7)

                        fitters[-1].SetFixRatio2GausSigma(0.6)
                        # fitters[-1].SetFixRatio2GausSigma(0.1)
                        # fitters[-1].SetBoundGaussianSigma(0.0005, 0.001)
                    status = fitters[-1].MassFitter(False)

                    if status != 1:  # fit failed
                        hYields.SetBinContent(iKStarBin+1, 0)
                        hYields.SetBinError(iKStarBin+1, 0)
                        hPurity.SetBinError(iKStarBin+1, 0)
                        hPurity.SetBinError(iKStarBin+1, 0)
                        hChi2.SetBinContent(iKStarBin+1, 0)
                        cMasses[iKStarBin // nPadsPerCanvas].cd(iKStarBin % nPadsPerCanvas + 1)
                        hMass.DrawClone()
                        continue

                    # fMass = fitters[-1].GetMassFunc()
                    fSgn = fitters[-1].GetSignalFunc()
                    fBkg = fitters[-1].GetBackgroundRecalcFunc()

                    hChi2.SetBinContent(iKStarBin + 1, fitters[-1].GetReducedChiSquare())

                    sgn, sgnUnc = ctypes.c_double(), ctypes.c_double()
                    bkg, bkgUnc = ctypes.c_double(), ctypes.c_double()
                    fitters[-1].Signal(2, sgn, sgnUnc)
                    fitters[-1].Background(2, bkg, bkgUnc)

                    hYields.SetBinContent(iKStarBin + 1, 0 if fSgn == None else sgn)
                    hYields.SetBinError(iKStarBin + 1, 0 if fSgn == None else sgnUnc)
                    hPurity.SetBinContent(iKStarBin + 1, 0 if fSgn == None else sgn.value/(sgn.value + bkg.value))
                    hPurity.SetBinError(iKStarBin + 1, 0 if fSgn == None else np.sqrt(bkg.value**2 * sgnUnc.value**2 + sgn.value**2 * bkgUnc.value **2) / (sgn.value + bkg.value)**2)

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
                hYields.Write()
                hPurity.Write()
                hChi2.Write()

                hDistr[event] = inFile.Get(f'{comb}/{event}/sgn/hCharmMassVsKStar0').ProjectionX()
                hDistr[event].SetName(f'h{event}PurityRew')
                hDistr[event].Rebin(round(kStarBW/hDistr[event].GetXaxis().GetBinWidth(1)))
                for iKStarBin in range(nKStarBins):
                    hDistr[event].SetBinContent(iKStarBin, hDistr[event].GetBinContent(iKStarBin + 1) * hPurity.GetBinContent(iKStarBin+1))
                hDistr[event].Write()
            
            hCF = CorrelationFunction(se=hDistr['SE'], me=hDistr['ME'], norm=pair.norm_range).get_cf()
            hCF.Write('hCFPurityRew')
            # break
        break
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
