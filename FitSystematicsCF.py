'''
Script to launch the fit on the correlation functions
obtained with systematic variations.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 FitSystematics.py cfg.yml

'''

import os
import argparse
import yaml
import ctypes
import numpy as np

# from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D, TSpline3
from ROOT import TFile, gInterpreter, TH1, TH1D, TH2D

from fempy import logger as log
from fempy.utils.io import Load

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

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Load yaml file
with open(f'{cfg["fitcfg"]}', "r") as stream:
    try:
        cfgfit = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml fit configuration could not be loaded. Is it properly formatted?')

outFileFits = cfgfit['ofilename']
if cfgfit['suffix']:
    outFileFits += f'_{cfgfit["suffix"]}'
dirFit, filenameFit = os.path.split(cfgfit['ofilename'])
if not os.path.isdir(dirFit + '/systfits'):
    os.makedirs(dirFit + '/systfits')
outFileFits = os.path.join(dirFit + '/systfits', filenameFit) + '_SystVarX.root' 

# Define the output file
oFileName = cfg['ofilename']
if cfg['suffix']:
    oFileName += f'_{cfg["suffix"]}'
oFileName += '.root'

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

nFits = len(cfgfit['fitcfs'])
fitParsList = [[] for nFit in range(nFits)] 
freeFitParsNames = [[] for nFit in range(nFits)] 
fitNames = [None for nFit in range(nFits)]
hDifferencesEntries = {}
hOfficialCFsEntries = {}

officialCFsFile = TFile(cfgfit['infile'])
hOfficialCFs = []
for iFit, fitName in enumerate(cfgfit['fitcfs']):
    hOfficialCF = Load(officialCFsFile, cfgfit['fitcfs'][iFit]['cfpath'])
    hOfficialCFs.append(hOfficialCF) 
    hOfficialCFsEntries[f"{cfgfit['fitcfs'][iFit]['fitname']}"] = []
    for iBin in range(hOfficialCF.GetNbinsX()):
        hOfficialCFsEntries[f"{cfgfit['fitcfs'][iFit]['fitname']}"].append(hOfficialCF.GetBinContent(iBin+1)) 

if cfgfit.get('evaluatediff'):
    # pick official differences
    hOfficialDifferences = {}

    # load file with fit results
    officialFitFile = TFile(f"{cfgfit['ofilename']}_{cfgfit['suffix']}.root")
    for iKey, key in enumerate(officialFitFile.GetListOfKeys()):
        fitNames[iKey] = key.GetTitle()
        hOfficialDifferences[key.GetTitle()] = Load(officialFitFile, f'{key.GetName()}/hDifference')
        hDifferencesEntries[f"{key.GetTitle()}"] = {}
        nBinsDifference = hOfficialDifferences[key.GetTitle()].GetNbinsX()
        for iBin in range(nBinsDifference):
            hDifferencesEntries[f"{key.GetTitle()}"][f"bin{iBin+1}"] = []

# fit the CFs with systematic variations
for iSystVar, systIdx in enumerate(cfg['systvars']):
    os.system(f"python3 fempy/FitCF.py {cfg['fitcfg']} --systvar {systIdx}")

    iSystVarFile = TFile(outFileFits.replace('X', str(systIdx)))

    # loop over the folders in the fit file
    for iKey, key in enumerate(iSystVarFile.GetListOfKeys()):
        iKeyFitPars = []
        iKeyFitParNames = []

        # load histo with fit parameters and histo checking 
        # whether a fit function parameter is fixed or not
        histoPars = Load(iSystVarFile, f'{key.GetName()}/hFitPars')
        histoFitCfg = Load(iSystVarFile, f'{key.GetName()}/hFreeFixPars')

        # save free fit parameters to a list
        for iBin in range(histoPars.GetNbinsX()):
            if(histoFitCfg.GetBinContent(iBin+3) == 1): # the first 2 bins contain the
                                                        # lower and upper fit ranges
                iKeyFitPars.append(histoPars.GetBinContent(iBin+1))
                iKeyFitParNames.append(histoPars.GetXaxis().GetBinLabel(iBin+1))

        fitParsList[iKey].append(iKeyFitPars)
        freeFitParsNames[iKey] = iKeyFitParNames

        if cfgfit.get('evaluatediff'):
            hDifference = Load(iSystVarFile, f'{key.GetName()}/hDifference')
            for iBin in range(nBinsDifference):
                hDifferencesEntries[f"{key.GetTitle()}"][f"bin{iBin+1}"].append(hDifference.GetBinContent(iBin+1)) 

# sanity check: compare mean of differences obtained with 
# systematic cuts with difference obtained with official cuts
for iFitName in fitNames:
    hDifferenceComparisons = TH2D(f"hDifferenceComp_{iFitName}", f"hDifferenceComp_{iFitName}", 
                                  nBinsDifference, 
                                  hOfficialDifferences[iFitName].GetBinLowEdge(1), 
                                  hOfficialDifferences[iFitName].GetBinLowEdge(nBinsDifference) + 
                                  hOfficialDifferences[iFitName].GetBinWidth(nBinsDifference), 
                                  len(cfg['systvars']), 0, len(cfg['systvars']))
    hDifferenceComparisons.SetStats(0)
    hDifferenceComparisons.GetYaxis().SetLabelSize(50)
    for iSystVar, systVar in enumerate(cfg['systvars']):
        hDifferenceComparisons.GetYaxis().SetBinLabel(iSystVar+1, f"{systVar}")
        for iBin in range(nBinsDifference):
            hDifferenceComparisons.SetBinContent(iBin+1, iSystVar+1, 
                                                 hOfficialDifferences[iFitName].GetBinContent(iBin+1) -
                                                 hDifferencesEntries[f"{iFitName}"][f"bin{iBin+1}"][iSystVar])

    oFile.mkdir(f'{iFitName}')
    oFile.cd(f'{iFitName}')
    hDifferenceComparisons.Write()

# systematics of the difference
hDifferencesSystUnc = []
for iFitName, fitName in enumerate(fitNames):
    hDifferenceWithSyst = hOfficialDifferences[fitName].Clone()
    hSystUnc = hOfficialDifferences[fitName].Clone()
    for iBin in range(nBinsDifference):
        hDifferenceWithSyst.SetBinError(iBin+1, np.std(hDifferencesEntries[f"{fitName}"][f"bin{iBin+1}"]))
        hSystUnc.SetBinContent(iBin+1, hDifferenceWithSyst.GetBinError(iBin+1))

    oFile.cd(f'{fitName}')
    hOfficialDifferences[fitName].Write("hDifferenceStat")
    hDifferenceWithSyst.Write("hDifferenceSyst")
    hSystUnc.Write("hSystUnc")

freeFitParsDicts = []
for iFit in range(nFits):
    freeParsHisto = {}
    for iPar in range(len(fitParsList[iFit][0])):
        parVars = []
        for fitParsVariation in fitParsList[iFit]:
            parVars.append(fitParsVariation[iPar])

        freeParsHisto[f'{freeFitParsNames[iFit][iPar]}'] = parVars
    freeFitParsDicts.append(freeParsHisto)

for iFitName in range(len(fitNames)):
    oFile.cd(f'{fitNames[iFitName]}')
    for iFitPar in freeFitParsDicts[iFitName]:
        histoPar = TH1D(f'h_{iFitPar}', f'h_{iFitPar}', 1000, 0.9 * min(freeFitParsDicts[iFitName][iFitPar]), 
                                                              1.1 * max(freeFitParsDicts[iFitName][iFitPar]))
        for iParValue in freeFitParsDicts[iFitName][iFitPar]:
            histoPar.Fill(iParValue)
        histoPar.Write()

oFile.Close()
print(f'\noutput saved in {oFileName}')
