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

# from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D, TSpline3
from ROOT import TFile, gInterpreter, TH1, TH1D

from fempy import logger as log
from fempy.utils.io import Load

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--systvar', default='0')
parser.add_argument('--debug', default=False, action='store_true')
parser.add_argument('--debugfit', default=False, action='store_true')
parser.add_argument('--debugdraw', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)
if args.debugfit:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_FIT 1')
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

excludeSystVars = cfg.get('exclsystvars', [])

nFits = len(cfgfit['fitcfs'])
fitParsList = [[] for nFit in range(nFits)] 
freeFitParsNames = [[] for nFit in range(nFits)] 
# hDifferenceList = [[None for iBin in range(0, cfgfit['fitcfs'][nFit][-1][1] / hOfficialCF)] for nFit in range(nFits)] 
fitNames = [None for nFit in range(nFits)]
# hOfficialCFs = [Load(TFile(cfgfit['infile']), cfgfit['fitcfs'][nfit]['cfpath']) for nFit in nFits]

print('Number of fits: ' + str(nFits))

# fit the CFs with systematic variations
for iSystVar in range(1, cfg['maxsystvar']+1):
    if iSystVar not in excludeSystVars:
        os.system(f"python3 fempy/FitCF.py {cfg['fitcfg']} --systvar {iSystVar}")

        iSystVarFile = TFile(outFileFits.replace('X', str(iSystVar)))
        # loop over the folders in the fit file
        for iKey, key in enumerate(iSystVarFile.GetListOfKeys()):
            print('KEYYYYYYYY ' + key.GetTitle())
            fitNames[iKey] = key.GetTitle()
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

            # try: # Check if there are ancestor histograms
            #     # Only load SE, the ME is the same for all
            #     hDifference = Load(inFile, f'{key.GetName()}/hDifference')
            #     hDifferenceList[iKey].append(hDifference)
            # except NameError: # Ancestors are not available, just move on
            #     log.warning(f"hDifference not found for syst var {iSystVar}")
            #     hDifferenceList.append(None)

freeFitParsDicts = []
for iFit in range(nFits):
    freeParsHisto = {}
    for iPar in range(len(fitParsList[iFit][0])):
        parVars = []
        for fitParsVariation in fitParsList[iFit]:
            parVars.append(fitParsVariation[iPar])

        freeParsHisto[f'{freeFitParsNames[iFit][iPar]}'] = parVars
    print("\n")
    print(freeParsHisto)
    print("\n")
    freeFitParsDicts.append(freeParsHisto)

for iFitName in range(len(fitNames)): 
    print(fitNames[iFitName])
    oFile.mkdir(f'{fitNames[iFitName]}')
    oFile.cd(f'{fitNames[iFitName]}')
    print('\n')
    print('\n')
    print('Set of fit parameters')
    print('\n')
    for iFitPar in freeFitParsDicts[iFitName]:
        histoPar = TH1D(f'h_{iFitPar}', f'h_{iFitPar}', 1000, 0.9 * min(freeFitParsDicts[iFitName][iFitPar]), 
                                                              1.1 * max(freeFitParsDicts[iFitName][iFitPar]))
        print('Low edge, min: ' + str(0.9 * min(freeFitParsDicts[iFitName][iFitPar])) + ' ' + 
                                  str(min(freeFitParsDicts[iFitName][iFitPar])) )
        print('Upp edge, max: ' + str(1.1 * max(freeFitParsDicts[iFitName][iFitPar])) + ' ' +
                                  str(max(freeFitParsDicts[iFitName][iFitPar]))   )
        for iParValue in freeFitParsDicts[iFitName][iFitPar]:
            histoPar.Fill(iParValue)
        histoPar.Write()

for iFitName in range(len(fitNames)): 
    oFile.mkdir(f'{fitNames[iFitName]}/differences')
    oFile.cd(f'{fitNames[iFitName]}/differences')

    # nBinsFitCF = round()

    # for iSystVar in range(1, cfg['maxsystvar']+1):
    #     if hDifferenceList[iFitName][iSystVar] is not None:


    #     double binWidth = this->fFitHist->GetBinWidth(1);
    #     int nBinsFitCF = static_cast<int>(this->fFitRangeMax/binWidth);
    # for

oFile.Close()
print(f'output saved in {oFileName}')
