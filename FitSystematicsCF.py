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

for iSystVar in range(1, cfg['maxsystvar']+1):
    if iSystVar not in excludeSystVars:
        os.system(f"python3 fempy/FitCF.py {cfg['fitcfg']} --systvar {iSystVar}")

oFile.cd()
fitParVars = []
for iSystVar in range(1, cfg['maxsystvar']+1):
    if iSystVar not in excludeSystVars:
        iVarFitPars = []
        iSystVarFile = TFile(cfg['outfitfiles'].replace('X', str(iSystVar)))
        for iKey in iSystVarFile.GetListOfKeys():
            iKeyFitPars = []
            histoPars = Load(iSystVarFile, f'{iKey.GetName()}/hFitPars')
            histoFitCfg = Load(iSystVarFile, f'{iKey.GetName()}/hFreeFixPars')
            for iBin in range(histoPars.GetNbinsX()):
                if(histoFitCfg.GetBinContent(iBin+3) == 1):
                    iKeyFitPars.append(histoPars.GetBinContent(iBin+1))
            iVarFitPars.append(iKeyFitPars)
            # histo.Write(f'hFitPars_{iKey.GetName()}_var{iSystVar}')
            nFreePars = len(iKeyFitPars)
        fitParVars.append(iVarFitPars)

fitParsVars = []
for iFreePar in range(nFreePars):
    for iFreePar


for iFitParVar in fitParVars:
    print(iFitParVar)
    print('\n')
fitPars = [[fitParVars[j][i] for j in range(len(fitParVars))] for i in range(len(fitParVars[0]))]
for iFitPar in fitPars:
    print(iFitPar)
    print('\n')


        
        # folder = Load(iSystVarFile, 'lpiplus')
        # folder.Write()

oFile.Close()
print(f'output saved in {oFileName}')
