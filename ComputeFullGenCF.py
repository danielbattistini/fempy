import argparse
import yaml
import itertools
import sys
import numpy as np

from alive_progress import alive_bar

from ROOT import TFile, TGraphErrors, TF1, gInterpreter, TDatabasePDG, gROOT, TNtuple, gRandom, TH1D, TH2D, TCanvas, kOrange, TLatex
import fempy
from fempy import CorrelationFunction

gInterpreter.ProcessLine('#include "fempy/utils/functions.h"')
gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter, GeneralCoulombLednicky


def Average(hist, xmin, xmax):
    firstBin = hist.GetXaxis().FindBin(xmin * 1.0001)
    lastBin = hist.GetXaxis().FindBin(xmax * 0.9999)

    avg = 0
    for iBin in range(firstBin, lastBin +1):
        avg += hist.GetBinContent(iBin) * hist.GetBinCenter(iBin)
    return avg/hist.Integral(firstBin, lastBin)


def StdDev(hist, xmin, xmax):
    firstBin = hist.GetXaxis().FindBin(xmin * 1.0001)
    lastBin = hist.GetXaxis().FindBin(xmax * 0.9999)
    mu = Average(hist, xmin, xmax)

    stdDev = 0
    for iBin in range(firstBin, lastBin +1):
        stdDev += hist.GetBinContent(iBin) * (hist.GetBinCenter(iBin) - mu)**2
    return (stdDev/hist.Integral(firstBin, lastBin))**0.5


def ApplyCenterOfGravity(hist, graph):
    nBins = hist.GetNbinsX()
    if nBins != graph.GetN():
        fempy.error("different number of bins")

    gCentered = graph.Clone(f'{hist.GetName()}_grav')
    for iBin in range(hist.GetNbinsX()):
        gCentered.SetPoint(iBin, graph.GetPointX(iBin), hist.GetBinContent(iBin+1))
        gCentered.SetPointError(iBin, graph.GetErrorX(iBin), hist.GetBinError(iBin+1))

    return gCentered


def IsLamParMatValid(lam_par):
    total = 0
    for _, h_lam in lam_par.items():
        for _, l_lam in h_lam.items():
            total += l_lam
    return abs(total - 1) < 1e-6


def IsSummedLamParMatValid(lam_par):
    total = 0
    for _, value in lam_par.items():
        total += value
    return abs(total - 1) < 1e-6


def SumLamPar(lam_par, treamtments):
    lam_par_summed = {
        'gen': 0,
        'flat': 0,
    }
    for hk, h_lam in lam_par.items():
        for lk, l_lam in h_lam.items():
            lam_par_summed[treamtments[hk][lk]] += l_lam
    return lam_par_summed


def RedMass(m1, m2):
    return m1 * m2 / (m1 + m2)


def Bootstrap(hist):
    hBootstrapped = hist.Clone(f'{hist.GetName()}_bootstrap')
    for iBin in range(hist.GetNbinsX()+2):
        mu = hist.GetBinContent(iBin)
        sigma = hist.GetBinError(iBin)

        hBootstrapped.SetBinContent(iBin, gRandom.Gaus(mu, sigma))
        hBootstrapped.SetBinError(iBin, sigma)
    return hBootstrapped


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair')
    parser.add_argument('--fitMass', default=False, action='store_true')
    parser.add_argument('--syst', default=False, action='store_true')
    parser.add_argument('--bs', type=int, default=0)
    parser.add_argument('-b', default=False, action='store_true')
    args = parser.parse_args()

    gROOT.SetBatch(args.b)

    # initialize
    regions = ['sgn', 'sbr']
    norm = [1500, 2000]
    kStarBW = 50 #MeV
    nKStarBins = round(3000 / kStarBW)
    nSystVar = 20 if args.syst else 1
    if args.pair == 'DstarK':
        fitRanges = [[10, 450], [10, 400], [10, 500]]
        inFileDistrData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileDistrMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = '/home/daniel/an/DstarK/2_luuksel/GenCF_nopc_kStarBW50MeV_fromq_sbcorr.root'
        config = '/home/daniel/phsw/fempy/cfg_gencf_DstarK_50MeV.yml'

        mLF = TDatabasePDG.Instance().GetParticle(321).Mass()

        weights1 = [0.78, 0.80, 0.77]
        radii1 = [0.86, 0.95, 0.79]
        radii2 = [2.03, 2.22, 1.91]
    elif args.pair == 'DstarPi':
        fitRanges = [[10, 350], [10, 300], [10, 400]]
        # fitRanges = [[10, 800], [10, 600], [10, 1000]]
        inFileDistrData = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileDistrMC = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_mcgp_nopc_kStarBW50MeV_true.root')
        oFileName = '/home/daniel/an/DstarPi/20_luuksel/GenCF_nopc_kStarBW50MeV.root'
        config = '/home/daniel/phsw/fempy/cfg_gencf_DstarPi_50MeV.yml'
        mLF = TDatabasePDG.Instance().GetParticle(211).Mass()
        weights1 = [0.66, 0.69, 0.64]
        radii1 = [0.97, 1.06, 0.89]
        radii2 = [2.52, 2.88, 2.32]

    # nSystVar = 5 #! debug
    if not args.syst:
        fitRanges = [fitRanges[0]]
        weights1 = [weights1[0]]
        radii1 = [radii1[0]]
        radii2 = [radii2[0]]
    redMass = RedMass(TDatabasePDG.Instance().GetParticle(413).Mass(), mLF)
    oFile = TFile(oFileName, 'recreate')

    with open(config, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    # compute raw CF of data
    for comb in ['sc', 'oc']:
        oFile.mkdir(comb)
        oFile.cd(comb)

        fempy.info(f'Bootstrapping the {comb} pairs')
        
        # Compute MC CF
        hSEMC = inFileDistrMC.Get(f'{comb}/SE/sgn/hCharmMassVsKStar0').ProjectionX('hSEMC')
        hMEMC = inFileDistrMC.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hMEMC')
        rebinFactor = round(kStarBW/hSEMC.GetBinWidth(1))
        hSEMC.Rebin(rebinFactor)
        hMEMC.Rebin(rebinFactor)
        cfMCOrig = CorrelationFunction(se=hSEMC, me=hMEMC, norm=norm)

        # Compute Gravities of bins
        hME = inFileDistrData.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hME')
        xCF = [Average(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
        xCFUnc = [StdDev(hME, iBin*kStarBW, (iBin+1)*kStarBW) for iBin in range(nKStarBins)]
        yCF = [1 for _ in range(nKStarBins)]
        yCFUnc = [0 for _ in range(nKStarBins)]
        gGravities = TGraphErrors(1)
        for iBin, (x, y, xUnc, yUnc) in enumerate(zip(xCF, yCF, xCFUnc, yCFUnc)):
            gGravities.SetPoint(iBin, x, y)
            gGravities.SetPointError(iBin, xUnc, yUnc)

        # Compute Data CF
        cfData = {}
        for region in regions:
            cfData[region] = {}
            for syst in range(nSystVar):
                hSEData = inFileDistrData.Get(f'{comb}/SE/{region}/hCharmMassVsKStar{syst}').ProjectionX(f'hSEData_{region}_{syst}')
                hMEData = inFileDistrData.Get(f'{comb}/ME/{region}/hCharmMassVsKStar{syst}').ProjectionX(f'hMEData_{region}_{syst}')
                hSEData.Rebin(rebinFactor)
                hMEData.Rebin(rebinFactor)
                cfData[region][syst] = CorrelationFunction(se=hSEData, me=hMEData, norm=norm)

        lamParMatr = {}
        for heavyContrib in cfg['heavy']:
            heavyKey = list(heavyContrib.keys())[0]
            heavyPurity = heavyContrib[heavyKey]['purity']
            heavyFrac = heavyContrib[heavyKey]['frac']
            lamParMatr[heavyKey] = {}
            for lightContrib in cfg['light']:
                lightKey = list(lightContrib.keys())[0]
                lightPurity = lightContrib[lightKey]['purity']
                lightFrac = lightContrib[lightKey]['frac']

                lamParMatr[heavyKey][lightKey] = lightFrac * lightPurity * heavyFrac * heavyPurity
        lamParCentr = SumLamPar(lamParMatr, cfg['treatment'])
        print("Lambda parameters for the central variation:")
        print(lamParMatr)
        print(lamParCentr)
        # sys.exit()
        if not IsSummedLamParMatValid(lamParCentr):
            print("Lambda parameters don't sum to 1. Exit!")
            sys.exit()

        # Variations on lambda parameters
        lamParSharp = dict(lamParCentr) # Explicit copy with dict()
        lamParSharp['flat'] *= 1.2
        lamParSharp['gen'] = 1 - lamParSharp['flat']

        lamParFlat = dict(lamParCentr)
        lamParFlat['flat'] *= 0.8
        lamParFlat['gen'] = 1 - lamParFlat['flat']

        lamPars = [lamParCentr, lamParFlat, lamParSharp]
        if not args.syst:
            lamPars = [lamPars[0]]

        # Perform multi trial
        tTrials = TNtuple('tTrials', 'trials', 'a0:a0unc:status:chi2ndf:iVar:w1:r1:r2:fitMax:lFlat:lGen')
        lCF = [[] for _ in range(60)]
        
        

        
        if not args.fitMass:
            varProduct = itertools.product(range(nSystVar), fitRanges, zip(radii1, radii2, weights1), lamPars)
            nTotalSystVars = nSystVar * len(fitRanges) * len(radii1) * len(lamPars) * (args.bs+1)
            with alive_bar(nTotalSystVars, force_tty=True) as bar:
                for iVar, (syst, fitRange, (radius1, radius2, weight1), lamPars) in enumerate(varProduct):
                    # Perform bootstrap for each syst variation
                    for iIter in range(args.bs+1): # In addition to the central variation that is always run
                        # if np.random.uniform(0, 1) > float(args.syst)/nTotalSystVars:
                        #     continue
                        if iIter > 0:
                            cfSgn = CorrelationFunction(cf=Bootstrap(cfData['sgn'][syst].get_cf()))
                            cfSbr = CorrelationFunction(cf=Bootstrap(cfData['sbr'][syst].get_cf()))
                            cfMC = CorrelationFunction(cf=Bootstrap(cfMCOrig.get_cf()))
                        else:
                            cfSgn = cfData['sgn'][syst]
                            cfSbr = cfData['sbr'][syst]
                            cfMC = cfMCOrig
                        # Fit the baseline
                        print("pirupiru")
                        fBaseLine = TF1('fBaseLine', '[0]', 0, 3000)
                        cfSgn.get_cf() # dummy call to compute the CF
                        hPurityDstar = TH1D("hPurityDstar", "", 60, 0, 3000)
                        for iBin in range(hPurityDstar.GetNbinsX()):
                            hPurityDstar.SetBinContent(iBin+1, 0.77)
                            hPurityDstar.SetBinError(iBin+1, 0.01)
                        purityDstar = CorrelationFunction(cf=hPurityDstar)

                        cfBkg = cfSgn/cfMC.get_cf()
                        ApplyCenterOfGravity(cfBkg.get_cf(), gGravities).Fit(fBaseLine, 'Q', '', 300, 1000)
                        blNorm = fBaseLine.GetParameter(0)

                        # Compute Gen CF
                        cfGen = (cfSgn/cfMC.get_cf()/blNorm - lamPars['flat'])/lamPars['gen']
                        
                        # fill lists for syst var
                        for iBin in range(cfGen.get_cf().GetNbinsX()):
                            lCF[iBin].append(cfGen.get_cf().GetBinContent(iBin+1))

                        if iIter == 0 and iVar == 0:
                            hCFGen = cfGen.get_cf().Clone()
                        # Fit the CF
                        gCFGen = ApplyCenterOfGravity(cfGen.get_cf(), gGravities)
                        fWeightedLL = TF1(f"fWeightedLL", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
                        fWeightedLL.FixParameter(0, radius1)
                        fWeightedLL.FixParameter(1, radius2)
                        fWeightedLL.FixParameter(2, weight1)
                        fWeightedLL.SetParameter(3, 0.1)
                        fWeightedLL.SetParLimits(3, -1, 1)
                        fWeightedLL.FixParameter(4, 0.)
                        fWeightedLL.FixParameter(5, 0.)
                        fWeightedLL.FixParameter(6, 1 if comb == 'sc' else -1)
                        fWeightedLL.FixParameter(7, redMass*1000)
                        
                        status = gCFGen.Fit(fWeightedLL, "SMRQ+0").Status()
                        chi2ndf = fWeightedLL.GetChisquare()/fWeightedLL.GetNDF()
                        scattLen = fWeightedLL.GetParameter(3)
                        scattLenUnc = fWeightedLL.GetParError(3)

                        tTrials.Fill(scattLen, scattLenUnc, status, chi2ndf, syst, weight1, radius1, radius2, fitRange[1], lamPars['flat'], lamPars['gen'])
                        bar()
        
        elif args.fitMass:
            # Perform multi trial
            print(f"Multi-trial info for the {comb} pair:")
            print("Variations on:")
            print("    n. of syst vars:   ", nSystVar)
            print("    n. of fit ranges:  ", len(fitRanges))
            print("    n. of radii:       ", len(radii1))
            print("    n. of lam-par vars:", len(lamPars))
            print("    n. of purity vars: ", 3 if args.syst else 1)

            varProduct = itertools.product(range(nSystVar), fitRanges, zip(radii1, radii2, weights1), lamPars)
            nPurityVar = 3 if args.syst else 1
            
            oFile.mkdir(f"{comb}/fits")
            oFile.cd(f"{comb}/fits")
            
            with alive_bar(nPurityVar*nSystVar * len(fitRanges) * len(radii1) * len(lamPars) * (args.bs+1), force_tty=True) as bar:
                for iVar, (syst, fitRange, (radius1, radius2, weight1), lamPars) in enumerate(varProduct):
                    # Perform bootstrap for each syst variation
                    # Make kStar slices
                    def ComputePurity(hist, kStarBW, event, variation=""):
                        nBins = round(3000/kStarBW)
                        hPurity = TH1D(f'h{event}Purity_{iVar}_{variation}', '', nBins, 0, 3000)
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

                    # Purity variations for SE
                    hSEMassVsKStar2D = inFileDistrData.Get(f'{comb}/SE/hCharmMassVsKStar{syst}')
                    hPuritySE = [ComputePurity(hSEMassVsKStar2D, kStarBW, 'SE')]
                    if args.syst:
                        hPuritySE.append(ComputePurity(hSEMassVsKStar2D, kStarBW, 'SE', '+1'))
                        hPuritySE.append(ComputePurity(hSEMassVsKStar2D, kStarBW, 'SE', '-1'))

                    hSEUnrew = inFileDistrData.Get(f'{comb}/SE/sgn/hCharmMassVsKStar{syst}').ProjectionX().Rebin(50)
                    oSEReweighted = [hPurity * hSEUnrew for hPurity in hPuritySE]
                    
                    # Purity variations for ME
                    hMEMassVsKStar2D = inFileDistrData.Get(f'{comb}/ME/hCharmMassVsKStar{syst}')
                    hPurityME = [ComputePurity(hMEMassVsKStar2D, kStarBW, 'ME')]
                    if args.syst:
                        hPurityME.append(ComputePurity(hMEMassVsKStar2D, kStarBW, 'ME', '+1'))
                        hPurityME.append(ComputePurity(hMEMassVsKStar2D, kStarBW, 'ME', '-1'))
                    hMEUnrew = inFileDistrData.Get(f'{comb}/ME/sgn/hCharmMassVsKStar{syst}').ProjectionX().Rebin(50)
                    oMEReweighted = [hPurity * hMEUnrew for hPurity in hPurityME]

                    for iPurityVar, (oSERew, oMERew) in enumerate(zip(oSEReweighted, oMEReweighted)):
                        oCFSBless = CorrelationFunction(se=oSERew, me=oMERew, norm=norm)
                        for iIter in range(args.bs+1): # In addition to the central variation that is always run
                            if iIter > 0:
                                cfsgn=CorrelationFunction(cf=Bootstrap(oCFSBless.get_cf()))
                                cfMC = CorrelationFunction(cf=Bootstrap(cfMCOrig.get_cf()))
                            else:
                                cfSgn = oCFSBless
                                cfMC = cfMCOrig

                            # Fit the baseline
                            fBaseLine = TF1('fBaseLine', '[0]', 0, 3000)
                            cfSgn.get_cf() # dummy call to compute the CF
                            # print('---------', cfMC.get_cf().GetNbinsX())
                            cfBkg = cfSgn/cfMC.get_cf()

                            ApplyCenterOfGravity(cfBkg.get_cf(), gGravities).Fit(fBaseLine, 'Q', '', 300, 600)
                            blNorm = fBaseLine.GetParameter(0)
                            blNormUnc = fBaseLine.GetParError(0)
                            hBL = cfSgn.get_cf().Clone('hBL')
                            for iBin in range(hBL.GetNbinsX()):
                                hBL.SetBinContent(iBin+1, blNorm)
                                hBL.SetBinError(iBin+1, blNormUnc)


                            # Compute Gen CF
                            cfGen = (cfSgn/cfMC.get_cf()/blNorm - lamPars['flat'])/lamPars['gen']
                            hGenCF = cfGen.get_cf()

                            for iBin in range(hGenCF.GetNbinsX()):
                                cfMJVal = cfMC.get_cf().GetBinContent(iBin+1)
                                cfSgnVal = cfSgn.get_cf().GetBinContent(iBin+1)
                                cfUnc = (1./lamPars['gen']/cfMJVal/blNorm * cfSgn.get_cf().GetBinError(iBin+1)) ** 2
                                cfUnc += (cfSgnVal/lamPars['gen']/blNorm/cfMJVal/cfMJVal * cfMC.get_cf().GetBinError(iBin+1)) ** 2
                                cfUnc += (cfSgnVal/lamPars['gen']/blNorm/blNorm/cfMJVal * blNormUnc) ** 2
                                # print(f"class: {hGenCF.GetBinError(iBin+1):.3f}    standalone: {cfUnc** 0.5 : .3f}")
                                hGenCF.SetBinError(iBin+1, cfUnc**0.5)


                            # cfGen = (cfSgn/cfMC.get_cf()/hBL - lamPars['flat'])/lamPars['gen']
                            if iIter == 0 and iVar == 0:
                                hCFGen = cfGen.get_cf().Clone()
                            # Fit the CF
                            gCFGen = ApplyCenterOfGravity(cfGen.get_cf(), gGravities)
                            fWeightedLL = TF1(f"fWeightedLL", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
                            fWeightedLL.FixParameter(0, radius1)
                            fWeightedLL.FixParameter(1, radius2)
                            fWeightedLL.FixParameter(2, weight1)
                            fWeightedLL.SetParameter(3, 0.1)
                            fWeightedLL.SetParLimits(3, -1, 1)
                            fWeightedLL.FixParameter(4, 0.)
                            fWeightedLL.FixParameter(5, 0.)
                            fWeightedLL.FixParameter(6, 1 if comb == 'sc' else -1)
                            fWeightedLL.FixParameter(7, redMass*1000)
                            
                            status = gCFGen.Fit(fWeightedLL, "SMRQ+0").Status()
                            chi2ndf = fWeightedLL.GetChisquare()/fWeightedLL.GetNDF()
                            scattLen = fWeightedLL.GetParameter(3)
                            scattLenUnc = fWeightedLL.GetParError(3)

                            tTrials.Fill(scattLen, scattLenUnc, status, chi2ndf, syst, weight1, radius1, radius2, fitRange[1], lamPars['flat'], lamPars['gen'])

                            # Draw fit
                            cFit = TCanvas(f'cFit_var{iVar}_iter{iIter}_iPurityVar{iPurityVar}', '', 600, 600)
                            cFit.DrawFrame(0, 0.8, fitRange[1], 1.2, fempy.utils.TranslateToLatex(';__kStarMeV__;__C__'))
                            gCFGen.Draw("same")

                            fWeightedLL.Draw("same")
                            
                            fWeightedLLCoulombOnly = TF1(f"fWeightedLLCoulombOnly", WeightedCoulombLednicky, fitRange[0], fitRange[1], 8)
                            fWeightedLLCoulombOnly.FixParameter(0, radius1)
                            fWeightedLLCoulombOnly.FixParameter(1, radius2)
                            fWeightedLLCoulombOnly.FixParameter(2, weight1)
                            fWeightedLLCoulombOnly.FixParameter(3, 0)
                            fWeightedLLCoulombOnly.FixParameter(4, 0.)
                            fWeightedLLCoulombOnly.FixParameter(5, 0.)
                            fWeightedLLCoulombOnly.FixParameter(6, 1 if comb == 'sc' else -1)
                            fWeightedLLCoulombOnly.FixParameter(7, redMass*1000)
                            fWeightedLLCoulombOnly.SetLineColor(kOrange+4)
                            fWeightedLLCoulombOnly.SetLineStyle(9)
                            fWeightedLLCoulombOnly.Draw("same")

                            tl = TLatex()
                            tl.SetNDC()
                            tl.SetTextSize(0.035)
                            tl.SetTextFont(42)
                            step = 0.05
                            tl.DrawLatex(0.6, 0.85 - 0 * step, f'ivar: {iVar}')
                            tl.DrawLatex(0.6, 0.85 - 1 * step, f'a_{{0}} = {scattLen:.3f} #pm {scattLenUnc:.3f}')
                            tl.DrawLatex(0.6, 0.85 - 2 * step, f'#chi^{{2}}/ndf = {chi2ndf:.2f}')

                            cFit.Write()


                            bar()

        oFile.cd(comb)
        
        tTrials.Write()
        cfMCOrig.get_cf().Write('hCFMC')
        for region in regions:
            for syst in range(nSystVar):
                cfData[region][syst].get_cf().Write(f'hCFData_{region}_{syst}')
        gGravities.Write('gGravities')

        if len(lCF[0]) > 0:
            yCF = [np.average(lCF[iBin]) for iBin in range(60)]
            yCFUnc = [np.std(lCF[iBin]) for iBin in range(60)]
        else:
            yCF = [cfGen.get_cf().GetBinContent(iBin+1) for iBin in range(60)]
            yCFUnc = [cfGen.get_cf().GetBinError(iBin+1) for iBin in range(60)]

        gGenCFTot = TGraphErrors(1)
        for iPoint in range(60):
            hCFGen.SetBinContent(iBin+1, yCF[iBin])
            gGenCFTot.SetPoint(iPoint, gCFGen.GetPointX(iPoint), yCF[iPoint])
            gGenCFTot.SetPointError(iPoint, gCFGen.GetErrorX(iPoint), yCFUnc[iPoint])

        ApplyCenterOfGravity(hCFGen, gGravities).Write('gGenCF')
        gGenCFTot.Write('gGenCFTot')
        hCFGen.Write('hGenCF')
        break
    oFile.Close()
    print(f'Output saved in {oFileName}')
