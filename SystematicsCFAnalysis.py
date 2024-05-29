'''
Script to compute the systematic uncertainties on the CF
bin values.
The output file is SystUncCF_suffix.root

Usage:
python3 SystematicsCFAnalysis.py cfg.yml

'''
import os
import argparse
import yaml
import math
from itertools import chain

from ROOT import TFile, TF1, TH1D, TH2D

from fempy import logger as log
from fempy.utils.io import Load, GetKeyNames

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

# check if the CFs for systematic variations have already been calculated, 
# else compute them
oFileCFBaseName = 'RawSystCF' + f'_{cfg["suffix"]}'
oFileCFName = os.path.join(cfg['odir'], oFileCFBaseName + '.root')
if os.path.exists(os.path.expanduser(oFileCFName)):
    print('Systematics CFs already computed!')
    inFileCFSystVar = TFile(oFileCFName)
else:
    print('Systematics CFs not yet computed!')
    os.system(f"python3 ComputeSystematicsCF.py {args.cfg}")
    print('Systematics CFs computed!')
    inFileCFSystVar = TFile(oFileCFName)
combs = [key.GetName() for key in inFileCFSystVar.GetListOfKeys() if key.GetName()[0] == 'p']

# Define the output file
oFileBaseName = 'SystUncCF'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)
    
# Load the CF selected for the analysis
inFileCFSelected = TFile(cfg['infileCFselected'])

# TH2D to store deviations from selected CF
systVars = list(chain.from_iterable(cfg['suffixessyst']))
nSystVars =  max(systVars) - min(systVars) + 1
nBinsKStar = len(range(cfg['systevalmaxbin']))
uppEdgeKStar = cfg['systevalmaxbin']*cfg['binwidth']

# Compute the correlation functions for all systematic variations
for comb in combs:
    print(f'Picking pair {comb}')
    oFile.mkdir(f'{comb}')
    oFile.cd(f'{comb}')

    hResiduals = TH2D("hResiduals", "Syst vars residuals", nBinsKStar, 0, uppEdgeKStar, 
                      nSystVars, min(systVars), max(systVars))
    hResiduals.SetStats(0)

    # list of histograms each containing the CF entries for
    # all variations for a specific bin 
    hCFSelected = Load(inFileCFSelected, f'{comb}/sgn/hCFrew')
    hCFBinEntries = []
    for iBin in range(cfg['systevalmaxbin']):
        hCFBinEntries.append(TH1D(f"h{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV", 
                                  f"h{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV", 
                                  cfg['systhistobins'], hCFSelected.GetBinContent(iBin+1)-cfg['systhistointerval'], 
                                  hCFSelected.GetBinContent(iBin+1)+cfg['systhistointerval']))
    
    for iSystVar in systVars:
        print(f'Picking syst variation number {iSystVar}')
        hSystVar = Load(inFileCFSystVar, f'{comb}/var{iSystVar}/hCFrew_{iSystVar}')
        print(f'Picked syst variation number {iSystVar}')

        hResiduals.GetYaxis().SetBinLabel(iSystVar-min(systVars)+1, str(iSystVar))
        for iBin in range(nBinsKStar):
            hCFBinEntries[iBin].Fill(hSystVar.GetBinContent(iBin+1))
            hResiduals.Fill(cfg['binwidth']*(iBin+1)+cfg['binwidth']/2, iSystVar,
                            (hCFSelected.GetBinContent(iBin+1)-hSystVar.GetBinContent(iBin+1)) / 
                            hCFSelected.GetBinContent(iBin+1))
 
    hCFYields = TH1D("hCFYields", "hCFYields", nBinsKStar, 0, uppEdgeKStar)
    hCFMeans = TH1D("hCFmeans", "hCFmeans", nBinsKStar, 0, uppEdgeKStar)
    hSystUnc = TH1D("hSystUnc", "hSystUnc", nBinsKStar, 0, uppEdgeKStar)
    hRelSystUnc = TH1D("hRelSystUnc", "hRelSystUnc", nBinsKStar, 0, uppEdgeKStar)
    hRelStatUnc = TH1D("hRelStatUnc", "hRelStatUnc", nBinsKStar, 0, uppEdgeKStar)
    hRatioStatSystUnc = TH1D("hRatioStatSystUnc", "hRatioStatSystUnc", nBinsKStar, 0, uppEdgeKStar)
    hCFWithStatSystUnc = TH1D("hCFWithStatSystUnc", "hCFWithStatSystUnc", nBinsKStar, 0, uppEdgeKStar)
    
    print('Evaluating the single bins ...')
    oFile.mkdir(f'{comb}/binsfits')
    for iBin, ihCFbin in enumerate(hCFBinEntries):
        oFile.cd(f'{comb}/binsfits')
        gaus = TF1(f"gaus_{iBin}", "gaus", ihCFbin.GetBinLowEdge(1), 
                   ihCFbin.GetBinLowEdge(ihCFbin.GetNbinsX()) + ihCFbin.GetBinWidth(1))
        gaus.SetParameter(1, hCFSelected.GetBinContent(iBin+1))
        ihCFbin.Fit(gaus, "SMRL+", "")
        ihCFbin.Write(f"hCFbin{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV")
        hCFYields.SetBinContent(iBin+1, gaus.GetParameter(0))
        hCFMeans.SetBinContent(iBin+1, gaus.GetParameter(1))
        hSystUnc.SetBinContent(iBin+1, gaus.GetParameter(2))
        hRelSystUnc.SetBinContent(iBin+1, hSystUnc.GetBinContent(iBin+1)/hCFSelected.GetBinContent(iBin+1))
        hRelStatUnc.SetBinContent(iBin+1, hCFSelected.GetBinError(iBin+1)/hCFSelected.GetBinContent(iBin+1))
        hRatioStatSystUnc.SetBinContent(iBin+1, hCFSelected.GetBinError(iBin+1)/hSystUnc.GetBinContent(iBin+1))
        hCFWithStatSystUnc.SetBinContent(iBin+1, hCFSelected.GetBinContent(iBin+1))
        hCFWithStatSystUnc.SetBinError(iBin+1, math.sqrt(hCFSelected.GetBinError(iBin+1)**2 + 
                                                       hSystUnc.GetBinContent(iBin+1)**2 ) )
    
    print('Writing histos with error informations ...')
    oFile.mkdir(f'{comb}/unc')
    oFile.cd(f'{comb}/unc')
    hCFYields.Write()
    hCFMeans.Write()
    hSystUnc.Write()
    hRelSystUnc.Write()
    hRelStatUnc.Write()
    hRatioStatSystUnc.Write()
    hCFWithStatSystUnc.Write()
    hResiduals.Write()

oFile.Close()
print(f'output saved in {oFileName}')
