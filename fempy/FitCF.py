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

from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D, TSpline3

from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import ChangeUnits

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--systvar', default='0')
parser.add_argument('--debug', default=False, action='store_true')
parser.add_argument('--debugfit', default=False, action='store_true')
parser.add_argument('--debugcombfit', default=False, action='store_true')
parser.add_argument('--debugdraw', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)
if args.debugfit:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_FIT 1')
if args.debugcombfit:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_COMBFIT 1')
if args.debugdraw:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_DRAW 1')

gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CombinedFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/DrawFitFuncts.hxx"')
from ROOT import CorrelationFitter, CombinedFitter, DrawFitFuncts

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Load input file with data CF
inFile = TFile(cfg['infilesyst'] if args.systvar != '0' else cfg['infile'])

# Define the output file
oFileName = cfg['ofilename']
if cfg['suffix']:
    oFileName += f'_{cfg["suffix"]}'
if args.systvar != '0':
    dir, filename = os.path.split(cfg['ofilename'])
    if not os.path.isdir(dir + '/systfits'):
        os.makedirs(dir + '/systfits')
    oFileName = os.path.join(dir + '/systfits', filename) + '_SystVar' + args.systvar 
oFileName += '.root'

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

fitters = []
drawFits = []
modelsBaselineIdxs = []
combPars = [[] for iFit in range(len(cfg['fitcfs']))]
modelsOnBaseline = [[] for iFit in range(len(cfg['fitcfs']))]
modelsShifts = [[] for iFit in range(len(cfg['fitcfs']))]
modelsMultNorm = [[] for iFit in range(len(cfg['fitcfs']))]
modelsMultGlobNorm = [[] for iFit in range(len(cfg['fitcfs']))]
modelsSaveSubComps = [[] for iFit in range(len(cfg['fitcfs']))]
modelsNormsSubComps = [[] for iFit in range(len(cfg['fitcfs']))]
modelsNormsSubCompsLabels = [[] for iFit in range(len(cfg['fitcfs']))]
modelsSubCompsMothers = [[] for iFit in range(len(cfg['fitcfs']))]
modelsSubComps = [[] for iFit in range(len(cfg['fitcfs']))]
modelsLegLabels = [[] for iFit in range(len(cfg['fitcfs']))]

# for loop over the correlation functions
for iFit, fitcf in enumerate(cfg['fitcfs']):

    # change unity of measure of histograms from GeV to MeV
    fitHisto = ChangeUnits(Load(inFile, fitcf['cfpath'] if args.systvar == '0'
                           else fitcf['cfsystpath'].replace('X', args.systvar)), 1000)

    # fit range
    fitters.append(CorrelationFitter(fitHisto, fitcf['fitrange']))
        
    # drawing class constructor
    drawRange = fitcf.get('drawrange', fitcf['fitrange'])
    drawFits.append(DrawFitFuncts(fitHisto, drawRange[0], drawRange[1]))

    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    baselineIdx = -1
    nPreviousPars = 0

    modelsLegLabels[iFit].append(fitcf['datalabel'])
    modelsLegLabels[iFit].append(fitcf['fitfunclabel'])
    
    # for loop over the functions entering in the model
    for iTerm, term in enumerate(fitcf['model']):
        modelsLegLabels[iFit].append(term.get('legentry', ''))
        modelsOnBaseline[iFit].append(term.get('onbaseline', 0))
        modelsShifts[iFit].append(term.get('shift', 0))
        modelsMultNorm[iFit].append(term.get('multnorm', 0))
        modelsMultGlobNorm[iFit].append(term.get('multglobnorm', 0))

        if cfg.get('combined'):
            for iKey, key in enumerate(term['params']):
                if 'shared' in key:
                    combPars[iFit].append(iKey + nPreviousPars + 1)
        nPreviousPars = nPreviousPars + len(term['params'])

        if term.get('subcomps'):
            modelsSaveSubComps[iFit].append(iTerm)
            modelsNormsSubComps[iFit].append(term['normssubcomps'])
            modelsNormsSubCompsLabels[iFit].append(term['normssubcompslabels'])
            modelsSubCompsMothers[iFit].append(term['func'])
            modelsSubComps[iFit].append(term['subcomps'])
            for iSubComp in range(len(term['subcomps'])):
                print('Reading subcomps!')
                modelsOnBaseline.append(term['sub_onbaseline'][iSubComp])                
                modelsLegLabels.append(term['sub_legentry'][iSubComp])
                modelsShifts.append(term.get('sub_shifts', [0.]*len(term['subcomps']))[iSubComp])
                modelsMultNorm.append(term.get('sub_multnorm', [1]*len(term['subcomps']))[iSubComp])
                modelsMultGlobNorm.append(term.get('sub_multglobnorm', [1]*len(term['subcomps']))[iSubComp])

        if term.get('isbaseline'):
            drawFits[-1].SetBasIdx(iTerm, term['addmode'] == "*")
            baselineIdx = iTerm
            modelsBaselineIdxs.append(baselineIdx)
                
        if term.get('template'):
            drawFits[-1].AddFitCompName(term['template'])
            templFile = TFile(term['templfile'])
            splinedTempl = Load(templFile, term['templpath'])
            if isinstance(splinedTempl, TH1):
                splinedTempl = ChangeUnits(splinedTempl, 1000)
                if term.get('rebin'):
                    splinedTempl.Rebin(term['rebin'])
            initPars = [(name, *vals) for name, vals in term['params'].items()]  
            fitters[-1].Add(term['template'], splinedTempl, initPars, term['addmode'], term.get('relweightcomp', 0))
            cSplinedTempl = TCanvas(f'c{term["template"]}', '', 600, 600)
            fitters[-1].DrawSpline(cSplinedTempl, splinedTempl)
            drawFits[-1].AddSplineHisto(splinedTempl)
            oFile.cd(fitcf['fitname'])
            cSplinedTempl.Write()
        
        elif term.get('func'):
            drawFits[-1].AddFitCompName(term['func'])
            if term.get('subcomps'):
                for iSubComp in range(len(term['subcomps'])):
                    drawFits[-1].AddFitCompName(term['subcomps'][iSubComp])

            if term.get('fixparsfromfuncts'):
                histoFuncFiles = []
                histoFuncHistoPars = []
                for histoFuncFile, histoFuncHistoParsPath in zip(term['histofuncfile'], term['histofuncpath']):
                    histoFuncFiles.append(TFile(histoFuncFile))
                    histoFuncHistoPars.append(Load(histoFuncFiles[-1], histoFuncHistoParsPath))
                    oFile.cd(fitcf['fitname'])

            initPars = []
            for iKey, key in enumerate(term['params']):
                if "fromfunct" == term['params'][key][0]:
                    initPars.append((key, histoFuncHistoPars[term['params'][key][1]].GetBinContent(term['params'][key][2]), 0, -1))
                else:
                    initPars.append((key, term['params'][key][0], term['params'][key][1], 
                                     term['params'][key][2]))

            fitters[-1].Add(term['func'], initPars, term['addmode'], term.get('relweightcomp', 0))
    
    # perform the fit and save the result
    if fitcf.get('globnorm'):
        fitters[-1].AddGlobNorm(f'{fitcf["fitname"]}_globnorm', fitcf['globnorm'][0], fitcf['globnorm'][1], fitcf['globnorm'][2])    
        drawFits[-1].SetGlobNorm(True)    
    
    fitters[-1].BuildFitFunction()

if cfg.get('combined'):
    combFitter = CombinedFitter(fitters, combPars)
    combFitter.CombinedFit()
else: 
    for iModel in range(len(fitters)):
        fitters[iModel].Fit()

        if cfg['fitcfs'][iModel].get('bootstraptries'):
            oFile.mkdir(f"{cfg['fitcfs'][iModel]['fitname']}/bootstrap")
            oFile.cd(f"{cfg['fitcfs'][iModel]['fitname']}/bootstrap")
            parBTDistros = fitters[-1].Bootstrap(cfg['fitcfs'][iModel]['bootstraptries'])
            for parBTDistro in parBTDistros:
                print('Integral of histo ' + parBTDistro.GetName() + ': ' + str(parBTDistro.Integral()))
                parBTDistro.Write()
            oFile.cd(cfg['fitcfs'][iModel]['fitname'])

        if cfg['fitcfs'][iModel].get('evaluatediff'):
            differenceHistos = fitters[-1].GetDifference(cfg['fitcfs'][iModel].get('bootstrapdifftries', 0))
            if len(differenceHistos)>1:
                oFile.mkdir(f"{cfg['fitcfs'][iModel]['fitname']}/bootstrap_diff")
                oFile.cd(f"{cfg['fitcfs'][iModel]['fitname']}/bootstrap_diff")
                for iDifferenceHisto in range(1, len(differenceHistos)):
                    differenceHistos[iDifferenceHisto].Write()
            oFile.cd(f"{cfg['fitcfs'][iModel]['fitname']}")
            differenceHistos[0].Write('hSubtraction')


for iModel in range(len(cfg['fitcfs'])):
    oFile.cd(cfg['fitcfs'][iModel]['fitname'])

    fitFunction = fitters[iModel].GetFitFunction()
    fitFunction.Write()
    fitHisto.Write()

    fitters[iModel].SaveFreeFixPars().Write()
    fitters[iModel].SaveFitPars().Write()
    fitters[iModel].PullDistribution().Write()

    hChi2DOF = TH1D('hChi2DOF', 'hChi2DOF', 1, 0, 1)
    if cfg.get('combined') is None:
        print('NOT COMBINED FIT')
        hChi2DOF.Fill(0.5, fitters[iModel].GetChi2Ndf())
        print('Chi2/NDF: ' + str(fitters[iModel].GetChi2Ndf()))
        hChi2DOFManual = TH1D('hChi2DOFManual', 'hChi2DOFManual', 1, 0, 1)
        hChi2DOFManual.Fill(0.5, fitters[iModel].GetChi2NdfManual())
        hChi2DOF.Write()
        hChi2DOFManual.Write()

    hAllCompsParHisto = fitters[iModel].SaveFitParsSplitComponents(modelsSubCompsMothers[iModel], modelsSaveSubComps[iModel], 
                                                                   modelsNormsSubComps[iModel], modelsSubComps[iModel], 
                                                                   modelsNormsSubCompsLabels[iModel])
    hAllCompsParHisto.Write()    

    cFit = TCanvas('cFit', '', 600, 600)
    drawFits[iModel].SetTotalFitFunc(fitFunction)
    drawFits[iModel].SetParHist(hAllCompsParHisto)

    if cfg['fitcfs'][iModel].get('drawsumcomps'):
        drawFits[iModel].EvaluateToBeDrawnComponents(modelsOnBaseline[iModel], modelsMultNorm[iModel],
                                                     modelsMultGlobNorm[iModel], modelsShifts[iModel],
                                                     modelsBaselineIdxs[iModel], cfg['fitcfs'][iModel]['drawsumcomps'])
        legLabels.extend(cfg['fitcfs'][iModel]['sumcompslegends'])
    else:        
        drawFits[iModel].EvaluateToBeDrawnComponents(modelsOnBaseline[iModel], modelsMultNorm[iModel],
                                                     modelsMultGlobNorm[iModel], modelsShifts[iModel],
                                                     modelsBaselineIdxs[iModel])

    drawFits[iModel].Draw(cFit, modelsLegLabels[iModel], cfg['fitcfs'][iModel]['legcoords'], cfg['fitcfs'][iModel]['linethick'])
    cFit.Write()

    if cfg['fitcfs'][iModel].get('isfitcf'):
        fitters[iModel].GetGenuine().Write("fGenuine")
        fitters[iModel].SaveScatPars().Write()

oFile.Close()
print(f'output saved in {oFileName}')