'''
Script to compute the purities of correlated particles.
The output file is RawPurities_suffix.root

Usage:
python3 ComputePurities.py cfg.yml

'''
import os
import argparse
import yaml
import numpy as np

from ROOT import TFile, TH1F, gInterpreter, TH1D, TCanvas, gPad
gInterpreter.ProcessLine('#include "../../../phsw/fempy/fempy/MassFitter.hxx"')
from ROOT import MassFitter
from fempy import logger as log
from fempy.utils.io import Load
import sys
sys.path.append('../../')
import fempy
from fempy.utils.__init__ import DivideCanvas

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--debug', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Define the output file
oFileBaseName = 'Purities'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

runSuffix = cfg['runsuffix']

chargedpart = {
    'p0': ['Particle0',''],
    'p1': ['Particle1','Anti']
}

neutralpart = {
    'p2': ['Particle2','v0'],
    'p3': ['Particle3','Antiv0']
}

# charged particles purity
if cfg['infilemc'] != '' and cfg['infilemc'] is not None:
    
    inFileMC = TFile(cfg['infilemc'])
    
    for cpart, fdcpart in chargedpart.items():
        
        # load MC histos relevant for purity
        identPartPtFile = f'HM{fdcpart[1]}TrkCutsMC{runSuffix}/HM{fdcpart[1]}TrkCutsMC{runSuffix}/IdentPartPt'
        hIdentPartPt = Load(inFileMC, identPartPtFile)
        corrPartPtFile = f'HM{fdcpart[1]}TrkCutsMC{runSuffix}/HM{fdcpart[1]}TrkCutsMC{runSuffix}/CorrParPt' 
        hCorrPartPt = Load(inFileMC, corrPartPtFile)
        
        # compute purity for each pT bin
        hCorrPartPt.Divide(hIdentPartPt)
        hPurity = hCorrPartPt
        hPurity.SetTitle(';#it{p}_{T} (GeV/#it{c});Purity')
        hPurity.SetName(f'hPurity_{cpart}')
        oFile.cd()
        hPurity.Write()
        
        # save purity with its error for each bin
        purities = []
        purityBinErrs = []
        for iBin in range(hPurity.GetNbinsX()):
            purities.append(hPurity.GetBinContent(iBin))
            nTotPart = hIdentPartPt.GetBinContent(iBin + 1)
            if(nTotPart != 0):
                binPurity = hPurity.GetBinContent(iBin + 1) 
                purityBinErr = np.sqrt( (binPurity * (1 - binPurity)) / nTotPart )
                hPurity.SetBinError(iBin + 1, purityBinErr )
                purityBinErrs.append(purityBinErr)
            else:
                purityBinErrs.append(0.)

        for npart, fdnpart in neutralpart.items():
            
            oFile.mkdir(f'{cpart}{npart[-1]}')
            oFile.cd(f'{cpart}{npart[-1]}')

            # load histo to compute averaged pT            
            kStarPtFile1 = f'HMResultsQA{runSuffix}/HMResultsQA{runSuffix}/PairQA/QA_{fdcpart[0]}_{fdnpart[0]}/'
            kStarPtFile2 = f'KstarPtSEPartOne_{fdcpart[0]}_{fdnpart[0]}'
            kStarPt = Load(inFileMC, kStarPtFile1 + kStarPtFile2)
            
            # compute yield of pions for each pT bin in k* range [0, 200] MeV/c2
            uppbin = kStarPt.GetXaxis().FindBin(0.2*0.9999) 
            hPt = kStarPt.ProjectionY('femto{fdcpart[0]}', 1, uppbin)
            hPt.Scale(1. / hPt.Integral())
            hPt.SetTitle(';#it{p}_{T} (GeV/#it{c});Counts')
            hPt.SetName(f'hPt_0_200')
            hPt.Write()
            
            # compute weights and their errors
            weights = []
            weightErrs = []
            if(hPurity.GetBinWidth(1) == hPt.GetBinWidth(1)):
                for iBin in range(hPurity.GetNbinsX()):
                    weights.append(hPt.GetBinContent(iBin))
                    weightErrs.append(hPt.GetBinError(iBin))
            else:
                for iBin in range(hPurity.GetNbinsX()):
                    
                    purityPt = hPurity.GetBinCenter(iBin)
                    weights.append(hPt.Interpolate(purityPt))
                    purityBin = hPt.FindBin(purityPt)
                    purityBinCenter = hPt.GetBinCenter(purityBin)
                    
                    if(purityBin <= purityBinCenter):
                        x0 = hPt.GetBinCenter(purityBin-1)
                        x1 = hPt.GetBinCenter(purityBin)
                        y0 = hPt.GetBinContent(purityBin-1)
                        y1 = hPt.GetBinContent(purityBin)
                        sy0 = hPt.GetBinError(purityBin-1)
                        sy1 = hPt.GetBinError(purityBin)
                        xFactor = (purityPt - x0) / (x1 - x0) 
                        sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                        weightErrs.append(sWeight)
                    
                    if(purityBin > purityBinCenter):
                        x0 = hPt.GetBinCenter(purityBin)
                        x1 = hPt.GetBinCenter(purityBin + 1)
                        y0 = hPt.GetBinContent(purityBin)
                        y1 = hPt.GetBinContent(purityBin + 1)
                        sy0 = hPt.GetBinError(purityBin)
                        sy1 = hPt.GetBinError(purityBin + 1)
                        xFactor = (purityPt - x0) / (x1 - x0) 
                        sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                        weightErrs.append(sWeight)            
            
            # compute purity
            avgPurity = 0
            for weight, purity in zip(weights, purities):
                avgPurity += weight * purity
            totWeights = sum(weights)
            avgPurity = avgPurity / totWeights
            
            # compute purity error
            gaussPurityErr = 0
            for weight, weighterr, purity, puritybinerr in zip(weights, weightErrs, purities, purityBinErrs):
                weiErr = purity * ( (1 / totWeights) + (weight / (totWeights**2) ) ) * weighterr
                purErr = weight * puritybinerr
                gaussPurityErr += weiErr**2 + purErr**2 
            avgPurityErr = np.sqrt(gaussPurityErr)        

            # save histo with average purity value for each k* bin
            SESampleFile1 = f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcpart[0]}_{fdnpart[0]}'
            SESampleFile2 = f'/SEDist_{fdcpart[0]}_{fdnpart[0]}'
            hSESample = Load(inFileMC, SESampleFile1 + SESampleFile2)
            hAvgPurity = hSESample.Clone('hAvgPurity')
            hAvgPurity.Reset('ICESM')
            hAvgPurity.SetTitle(';#it{k}^{*} (GeV/#it{c});Purity')
            for iBin in range(hAvgPurity.GetNbinsX()):
                hAvgPurity.SetBinContent(iBin, avgPurity) 
                hAvgPurity.SetBinError(iBin, avgPurityErr) 
            hAvgPurity.Write()

# neutral particles purity
if cfg['infiledata'] != '' and cfg['infiledata'] is not None:
    
    inFileData = TFile(cfg['infiledata'])
    nPartList = cfg['neutrpart']
    nPartNumber = range(len(nPartList))
    
    for nPart, number in zip(nPartList, nPartNumber):
        
        nPartType = nPart['type']
        nPartName = nPart['name']
        nSigma = nPart['nsigma']
        ptMins = nPart['ptMins']
        ptMaxs = nPart['ptMaxs']
        sgnFuncName = nPart['sgnfuncname']
        bkgFuncName = nPart['bkgfuncname']
        fitRange = nPart['fitrange']
        task = nPart['task']
        fitSettingList = nPart['fitsettings']
        cfgFilePath = nPart['cfgfile']

        oFile.mkdir(nPartType)

        # load TH2D with invariant mass spectra
        hInvMassPt = Load(inFileData, f'HM{nPartName}Cuts{runSuffix}/HM{nPartName}Cuts{runSuffix}/v0Cuts/InvMassPt')
        
        # create useful histos
        hInvMassWidthVsPt = TH1D('hInvMassWidthVsPt', ';#it{p}_{T} (GeV/#it{c});Width (GeV/#it{c})', len(ptMins), 
                                 np.array(ptMins + [ptMaxs[-1]], 'd'))
        hInvMassMeanVsPt = TH1D('hInvMassMeanVsPt', ';#it{p}_{T} (GeV/#it{c});Mean (GeV/#it{c})', len(ptMins), 
                                np.array(ptMins + [ptMaxs[-1]], 'd'))
        hPurityVsPt = TH1D('hPurityVsPt', ';#it{p}_{T} (GeV/#it{c});Purity', len(ptMins), 
                           np.array(ptMins + [ptMaxs[-1]], 'd'))
        hChi2VsPt = TH1D('hChi2VsPt', ';#it{p}_{T} (GeV/#it{c});Chi2/DOF', len(ptMins), 
                         np.array(ptMins + [ptMaxs[-1]], 'd'))
        
        # create list of fit configuration
        fitters = []
        cInvMass = TCanvas('cInvMass', '', 1200, 800)
        DivideCanvas(cInvMass, len(ptMins))
        
        for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
            
            cInvMass.cd(iPt+1)
            
            # project the TH2D to obtain the invariant mass spectra to fit 
            lowBin = hInvMassPt.GetXaxis().FindBin(ptMin*1.0001) 
            uppBin = hInvMassPt.GetXaxis().FindBin(ptMax*0.9999) 
            hMassProj = hInvMassPt.ProjectionY(f'hMass_{ptMin}_{ptMax}', lowBin, uppBin)
            hMassProj.SetTitle(f'{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}}')
            
            # execute fit with settings from config file
            fitters.append(MassFitter(hMassProj, sgnFuncName, bkgFuncName, fitRange[0], fitRange[1], cfgFilePath))
            fitters[-1].SetFitSettings(task, number, iPt)
            fitters[-1].Fit()
            fitters[-1].Draw(gPad, "data_minus_bkg")
                        
            # fill histos with fit information
            hInvMassWidthVsPt.SetBinContent(iPt+1, fitters[-1].GetWidth())
            hInvMassWidthVsPt.SetBinError(iPt+1, fitters[-1].GetWidthUnc())
            hInvMassMeanVsPt.SetBinContent(iPt+1, fitters[-1].GetMean())
            hInvMassMeanVsPt.SetBinError(iPt+1, fitters[-1].GetMeanUnc())
            hChi2VsPt.SetBinContent(iPt+1, fitters[-1].GetChi2Ndf())
            
            # compute purity
            sgnFunc = fitters[-1].GetSgnFunc()
            bkgFunc = fitters[-1].GetBkgFunc() 
            sgnMax = sgnFunc.GetMaximumX() 
            sgnWidth =  nSigma * fitters[-1].GetWidth() 
            uppIntEdge = sgnMax + sgnWidth
            lowIntEdge = sgnMax - sgnWidth
            sgnInt = sgnFunc.Integral(lowIntEdge, uppIntEdge)
            bkgInt = bkgFunc.Integral(lowIntEdge, uppIntEdge) 
            purity = sgnInt / (sgnInt + bkgInt)
            hPurityVsPt.SetBinContent(iPt+1, purity)
            
            # compute purity error
            sgnIntErr = fitters[-1].GetSignalUnc(nSigma, "sgn_int")  
            bkgIntErr = fitters[-1].GetBackgroundUnc(nSigma) 
            dSgn = 1 / (sgnInt + bkgInt) - sgnInt / ( (sgnInt + bkgInt)**2 )
            dBkg = sgnInt / ( (sgnInt + bkgInt)**2 )
            purityErr = np.sqrt( (dSgn * sgnIntErr)**2 + (dBkg * bkgIntErr)**2 )
            hPurityVsPt.SetBinError(iPt+1, purityErr)
            
            oFile.cd(nPartType)
            hMassProj.Write()

        hInvMassWidthVsPt.Write()
        hInvMassMeanVsPt.Write()
        hPurityVsPt.Write()
        hChi2VsPt.Write()
        cInvMass.Write()
        
oFile.Close()
print(f'output saved in {oFileName}')

