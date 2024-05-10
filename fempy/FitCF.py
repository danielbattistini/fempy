'''
Script to perform the fit on a correlation function.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 CorrelationFitter.py cfg.yml

'''

import os
import argparse
import yaml
import ctypes

from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/DrawFitFuncts.hxx"')
from ROOT import CorrelationFitter, DrawFitFuncts

from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import ChangeUnits

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

# Load input file with data and mc CF
inFileFit = TFile(cfg['infile'])

# Define the output file
oFileBaseName = cfg['ofilebasename']
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

modelFitters = []
drawFits = []

# for loop over the correlation functions
for iFit, fitcf in enumerate(cfg['fitcfs']):
    
    # change unity of measure of histograms from GeV to MeV
    fitHisto = ChangeUnits(Load(inFileFit, fitcf['cfpath']), 1000)
    if('rebin' in fitcf):
        fitHisto.Rebin(fitcf['rebin'])

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    if('rejectrange' in fitcf):
        lowRejectRange = fitcf['rejectrange'][0]
        uppRejectRange = fitcf['rejectrange'][1]
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange, lowRejectRange, uppRejectRange))
    else: 
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange))
        
    # drawing class constructor
    if('drawrange' in fitcf):
        lowDrawRange = fitcf['drawrange'][0]
        uppDrawRange = fitcf['drawrange'][1]
        drawFits.append(DrawFitFuncts(fitHisto, lowDrawRange, uppDrawRange))
    else: 
        drawFits.append(DrawFitFuncts(fitHisto, lowFitRange, uppFitRange))

    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    baselineIdx = -1
    linesThickness = fitcf['linethick']
    compsToFile = []
    onBaseline = []
    shifts = []
    multNorm = []
    multGlobNorm = []
    saveSubComps = []
    normsSubComps = []
    normsSubCompsLabels = []
    subComps = []
    legLabels = []
    legLabels.append(fitcf['datalabel'])
    legLabels.append(fitcf['fitfunclabel'])
    # for loop over the functions entering in the model
    for iTerm, term in enumerate(fitcf['model']):
        
        if('subcomps' in term):
            saveSubComps.append(iTerm-1)
            normsSubComps.append(term['normssubcomps'])
            normsSubCompsLabels.append(term['normssubcompslabels'])
            subComps.append(term['subcomps'])
            print('Reading subcomps!')
            for iSubComp in range(len(term['subcomps'])):
                onBaseline.append(term['sub_onbaseline'][iSubComp])                
                legLabels.append(term['sub_legentry'][iSubComp])
                shifts.append(term.get('sub_shifts', [0.]*len(term['subcomps']))[iSubComp])
                multNorm.append(term.get('sub_multnorm', [1]*len(term['subcomps']))[iSubComp])
                multGlobNorm.append(term.get('sub_multglobnorm', [1]*len(term['subcomps']))[iSubComp])
        
        onBaseline.append(term['onbaseline'])
        legLabels.append(term['legentry'])
        shifts.append(term.get('shift', 0.))
        multNorm.append(term.get('multnorm', 1))
        multGlobNorm.append(term.get('multglobnorm', 1))

        if('isbaseline' in term):
            drawFits[-1].SetBasIdx(iTerm)
            baselineIdx = iTerm
                
        if('template' in term):
            drawFits[-1].AddFitCompName(term['template'])
            templFile = TFile(term['templfile'])
            splinedTempl = Load(templFile, term['templpath'])
            drawFits[-1].AddSplineHisto(splinedTempl)
            if(isinstance(splinedTempl, TH1)):
                splinedTempl = ChangeUnits(splinedTempl, 1000)
                if('rebin' in term):
                    splinedTempl.Rebin(term['rebin'])
            initPars = [(key, term['params'][key][0], term['params'][key][1], 
                         term['params'][key][2]) for key in term['params']]    
            modelFitters[-1].Add(term['template'], splinedTempl, initPars, term['addmode'])
            cSplinedTempl = TCanvas(f'c{term["template"]}', '', 600, 600)
            modelFitters[-1].DrawSpline(cSplinedTempl, splinedTempl)
            oFile.cd(fitcf['fitname'])
            cSplinedTempl.Write()
        
        elif('spline' in term):
            drawFits[-1].AddFitCompName(term['spline'])
            histoFile = TFile(term['histofile'])
            toBeSplinedHisto = ChangeUnits(Load(histoFile, term['histopath']), term['changeunits'])
            drawFits[-1].AddSplineHisto(toBeSplinedHisto)
            initPars = []
            normPar = list(term['params'].keys())[0]
            initPars.append((normPar, term['params'][normPar][0], 
                             term['params'][normPar][1], term['params'][normPar][2]))
            for nKnot, xKnot in enumerate(term['params']['xknots']):
                initPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
            for nKnot, xKnot in enumerate(term['params']['xknots']):
                nBin = toBeSplinedHisto.FindBin(xKnot)
                yKnot = toBeSplinedHisto.GetBinContent(nBin)
                initPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*term['yboundperc'], yKnot + (yKnot/100)*term['yboundperc']])
            modelFitters[-1].Add(term['spline'], initPars, term['addmode'])
        
        elif('func' in term):
            drawFits[-1].AddFitCompName(term['func'])
            if('subcomps' in term):
                for iSubComp in range(len(term['subcomps'])):
                    drawFits[-1].AddFitCompName(term['subcomps'][iSubComp])

            if('fixparsfromfuncts' in term):
                histoFuncFiles = []
                histoFuncHistoPars = []
                for histoFuncFile, histoFuncHistoParsPath in zip(term['histofuncfile'], term['histofuncpath']):
                    histoFuncFiles.append(TFile(histoFuncFile))
                    histoFuncHistoPars.append(Load(histoFuncFiles[-1], histoFuncHistoParsPath))

            initPars = []
            for key in term['params']:
                if("fromfunct" == term['params'][key][0]):
                    initPars.append((key, histoFuncHistoPars[term['params'][key][1]].GetBinContent(term['params'][key][2]), 0, -1))
                else:
                    initPars.append((key, term['params'][key][0], term['params'][key][1], 
                                     term['params'][key][2]))

            modelFitters[-1].Add(term['func'], initPars, term['addmode'])
                
    # perform the fit and save the result
    if('globnorm' in fitcf):
        modelFitters[-1].AddGlobNorm('globnorm', fitcf['globnorm'][0], fitcf['globnorm'][1], fitcf['globnorm'][2])    
        drawFits[-1].SetGlobNorm(True)    
    
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    modelFitters[-1].Fit()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    fitFunction.Write()
    hChi2DOF = TH1D('hChi2DOF', 'hChi2DOF', 1, 0, 1)
    hChi2DOF.Fill(0.5, modelFitters[-1].GetChi2Ndf())
    print('Chi2 / DOF: ' + str(modelFitters[-1].GetChi2Ndf()))
    print('\n\n')
    hChi2DOFManual = TH1D('hChi2DOFManual', 'hChi2DOFManual', 1, 0, 1)
    hChi2DOFManual.Fill(0.5, modelFitters[-1].GetChi2NdfManual())
    hChi2DOF.Write()
    modelFitters[-1].SaveFreeFixPars().Write()
    modelFitters[-1].SaveFitPars().Write()
    modelFitters[-1].PullDistribution().Write()
    if('isfitcf' in fitcf):
        modelFitters[-1].SaveScatPars().Write()

    hAllCompsParHisto = modelFitters[-1].SaveFitParsSplitComponents(saveSubComps, normsSubComps, subComps, normsSubCompsLabels)
    hAllCompsParHisto.Write()    
    cFit = TCanvas('cFit', '', 600, 600)
    drawFits[-1].SetTotalFitFunc(fitFunction)
    drawFits[-1].SetParHist(hAllCompsParHisto)
    if('drawsumcomps' in fitcf):
        drawFits[-1].EvaluateToBeDrawnComponents(onBaseline, multNorm, multGlobNorm, shifts,
                                             baselineIdx, fitcf['drawsumcomps'])
    else:        
        drawFits[-1].EvaluateToBeDrawnComponents(onBaseline, multNorm, multGlobNorm, shifts,
                                             baselineIdx)
    
    drawFits[-1].Draw(cFit, legLabels, fitcf['legcoords'], linesThickness)
    
    if('isfitcf' in fitcf):
        modelFitters[-1].GetGenuine().Write("fGenuine")
    if('debug' in fitcf):
        modelFitters[-1].Debug()
        
    cFit.Write()

    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    modelFitters[-1].SaveFitPars().Write()
    modelFitters[-1].SaveFreeFixPars().Write()
    modelFitters[-1].PullDistribution().Write()
    if('isfitcf' in fitcf):
        modelFitters[-1].SaveScatPars().Write()
    fitFunction.Write()
    hChi2DOF.Write()
    hChi2DOFManual.Write()

    print('Getting components ...')
    for iTerm, term in enumerate(fitcf['model']):
        if('savetofile' in term):
            if(term['savetofile']):
                compsToFile.append([iTerm, term.get('onbaseline', 0), baselineIdx])

    print(compsToFile)
    print('Saved, now writing to file ...')
    compsToBeWritten = modelFitters[-1].GetComponents([idx[0] for idx in compsToFile], 
                                                      [onBas[1] for onBas in compsToFile],
                                                      baselineIdx)
    print(len(compsToBeWritten))
    for comp in compsToBeWritten:
        comp.Write("f" + comp.GetName())
    print('Saved components!')

oFile.Close()
print(f'output saved in {oFileName}')