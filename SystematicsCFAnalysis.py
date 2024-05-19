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

from ROOT import TFile, TF1, TH1D

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
inFileCFSystVar = TFile(cfg['infileCFsystvar'])

if 'combs' in cfg:
    combs = cfg['combs']
else:
    combs = ['p02', 'p03', 'p12', 'p13', 'p02_13', 'p03_12']

# Compute the correlation functions for all systematic variations
for comb in combs:
    print(f'Picking pair {comb}')
    oFile.mkdir(comb)
    oFile.mkdir(f'{comb}/binsfits')
    oFile.mkdir(f'{comb}/unc')
    oFile.cd(comb)

    # list of histograms each containing the CF entries for
    # all variations for a specific bin 
    hCFBinEntries = []
    for iBin in range(cfg['systevalnbins'][0], cfg['systevalnbins'][1]):
        hCFBinEntries.append(TH1D(f"h{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV", 
                                  f"h{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV", 
                                  cfg['systhistobins'], cfg['systhistolowedge'], 
                                  cfg['systhistouppedge']))
    
    for iSystVar in range(1, cfg['systnvariations']):
        print(f'Picking syst variation number {iSystVar}')
        hSystVar = Load(inFileCFSystVar, f'{comb}/CFsystvar/{iSystVar}/hCFrew_{iSystVar}')

        for iBin in range(cfg['systevalnbins'][0], cfg['systevalnbins'][1]):
            hCFBinEntries[iBin].Fill(hSystVar.GetBinContent(iBin+1))

    #for iBin, ihCFbin in enumerate(hCFBinEntries):
    #    oFile.cd(f'{comb}/binsfits')
    #    ihCFbin.Write(f"hCFbin{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV")
        
    nBins = cfg['systevalnbins'][1] - cfg['systevalnbins'][0] 
    lowEdgeHisto = cfg['systevalnbins'][0]*cfg['binwidth']
    uppEdgeHisto = cfg['systevalnbins'][1]*cfg['binwidth']
    hCFSelected = Load(inFileCFSelected, f'{comb}/sgn/hCFrew')
    hCFYields = TH1D("hCFYields", "hCFYields", nBins, lowEdgeHisto, uppEdgeHisto)
    hCFMeans = TH1D("hCFmeans", "hCFmeans", nBins, lowEdgeHisto, uppEdgeHisto)
    hSystUnc = TH1D("hSystUnc", "hSystUnc", nBins, lowEdgeHisto, uppEdgeHisto)
    hRelSystUnc = TH1D("hRelSystUnc", "hRelSystUnc", nBins, lowEdgeHisto, uppEdgeHisto)
    hRelStatUnc = TH1D("hRelStatUnc", "hRelStatUnc", nBins, lowEdgeHisto, uppEdgeHisto)
    hRatioStatSystUnc = TH1D("hRatioStatSystUnc", "hRatioStatSystUnc", nBins, lowEdgeHisto, uppEdgeHisto)
    hCFWithStatSystUnc = TH1D("hCFWithStatSystUnc", "hCFWithStatSystUnc", nBins, lowEdgeHisto, uppEdgeHisto)
    
    print('Evaluating the single bins ...')
    for iBin, ihCFbin in enumerate(hCFBinEntries):
        oFile.cd(f'{comb}/binsfits')
        gaus = TF1(f"gaus_{iBin}", "gaus", cfg['systhistolowedge'], cfg['systhistouppedge'])
        gaus.SetParameter(1, hCFSelected.GetBinContent(iBin+1))
        ihCFbin.Fit(gaus, "SMRL+0", "")
        ihCFbin.Write(f"hCFbin{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV")
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
    oFile.cd(f'{comb}/unc')
    hCFYields.Write()
    hCFMeans.Write()
    hSystUnc.Write()
    hRelSystUnc.Write()
    hRelStatUnc.Write()
    hRatioStatSystUnc.Write()
    hCFWithStatSystUnc.Write()

oFile.Close()
print(f'output saved in {oFileName}')
