'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import argparse
import yaml

from ROOT import TFile, TH1F

from fempy import logger as log
from fempy.utils.io import Load

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

regions = ['sgn']
combs = {
    'p02': 'Particle0_Particle2',
    'p03': 'Particle0_Particle3',
    'p12': 'Particle1_Particle2',
    'p13': 'Particle1_Particle3'
}

# Load input file with same- and mixed-event distributions
inFile = TFile(cfg['infile'])

# Define the output file
oFileBaseName = 'RawCF'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

hSE = {}
hME = {}
hMErew = {}
for comb, fdcomb in combs.items():
    oFile.mkdir(comb)
    oFile.cd(comb)
    
    hSE[comb] = {}
    hME[comb] = {}
    hMErew[comb] = {}
    for region in regions:
        pcsuffix = cfg['pcsuffix']
        hSE[comb][region] = Load(inFile, f'HMResults{pcsuffix}/HMResults{pcsuffix}/{fdcomb}/SEDist_{fdcomb}')
        hME[comb][region] = Load(inFile, f'HMResults{pcsuffix}/HMResults{pcsuffix}/{fdcomb}/MEDist_{fdcomb}')
        hSEmultk = Load(inFile, f'HMResults{pcsuffix}/HMResults{pcsuffix}/{fdcomb}/SEMultDist_{fdcomb}')
        hMEmultk = Load(inFile, f'HMResults{pcsuffix}/HMResults{pcsuffix}/{fdcomb}/MEMultDist_{fdcomb}')
        nbins = hMEmultk.ProjectionX().GetNbinsX()
        hMEreweightk = TH1F("MErewdistr", "MErewdistr", nbins, hMEmultk.GetXaxis().GetXmin(), hMEmultk.GetXaxis().GetXmax())
        for iBin in range(hMEmultk.ProjectionY().GetNbinsX()):
            hSEbinmult = hSEmultk.ProjectionX(f'{comb}SEdistr', iBin, iBin+1)
            hMEbinmult = hMEmultk.ProjectionX(f'{comb}MEdistr', iBin, iBin+1)
            hMEreweightk.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
        hMErew[comb][region] = hMEreweightk
   
# Sum pair and antipair
for comb in combs:
    hSE['p02_13'] = {}
    hSE['p03_12'] = {}
    hME['p02_13'] = {}
    hME['p03_12'] = {}
    hMErew['p02_13'] = {}
    hMErew['p03_12'] = {}

    for region in regions:
        hSE['p02_13'][region] = hSE['p02'][region] + hSE['p13'][region]
        hME['p02_13'][region] = hME['p02'][region] + hME['p13'][region]
        hMErew['p02_13'][region] = hMErew['p02'][region] + hMErew['p13'][region]
        
        hSE['p03_12'][region] = hSE['p03'][region] + hSE['p12'][region]
        hME['p03_12'][region] = hME['p03'][region] + hME['p12'][region]
        hMErew['p03_12'][region] = hMErew['p03'][region] + hMErew['p12'][region]

# Compute the CF and write to file
for comb in list(combs.keys()) + ['p02_13', 'p03_12']:
    for region in regions:
        rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
        print(rebin)
        hSE[comb][region].Rebin(rebin)
        hME[comb][region].Rebin(rebin)
        hMErew[comb][region].Rebin(rebin)

        if cfg['norm'] is None:  # if not specified, normalize to the yields
            norm = hME[comb][region].GetEntries() / hSE[comb][region].GetEntries()
            normrew = hMErew[comb][region].GetEntries() / hSE[comb][region].GetEntries()
        else:
            log.critical('Normalization method not implemented')

        hSE[comb][region].Sumw2()
        
        hCF = norm * hSE[comb][region] / hME[comb][region]
        hCFrew = normrew * hSE[comb][region] / hMErew[comb][region]

        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSE[comb][region].SetName('hSE')
        hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hME[comb][region].SetName('hME')
        hME[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        
        hMErew[comb][region].SetName('hMErew')
        hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hCF.SetName('hCF')
        hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCF.Write()
        
        hCFrew.SetName('hCFrew')
        hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCFrew.Write()

oFile.Close()
print(f'output saved in {oFileName}')
