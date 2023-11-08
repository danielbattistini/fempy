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
    cfg['pairs'][0]: 'Particle0_Particle2',
    cfg['pairs'][1]: 'Particle0_Particle3',
    cfg['pairs'][2]: 'Particle1_Particle2',
    cfg['pairs'][3]: 'Particle1_Particle3'
}

firstpair = cfg['pairs'][0] + cfg['pairs'][3]
secondpair = cfg['pairs'][1] + cfg['pairs'][2]

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
        hSE[comb][region] = Load(inFile, f'HMResults1001/HMResults1001;1/{fdcomb}/SEDist_{fdcomb}')
        hME[comb][region] = Load(inFile, f'HMResults1001/HMResults1001;1/{fdcomb}/MEDist_{fdcomb}')
        binwidth = cfg['multbinwidth']
        hSEmultk = Load(inFile, f'HMResults1001/HMResults1001;1/{fdcomb}/SEMultDist_{fdcomb}')
        hMEmultk = Load(inFile, f'HMResults1001/HMResults1001;1/{fdcomb}/MEMultDist_{fdcomb}')
        nbins = hMEmultk.ProjectionX().GetNbinsX()
        hMEreweightk = TH1F("MErewdistr", "MErewdistr", nbins, hMEmultk.ProjectionX().GetBinLowEdge(1), 
                            hMEmultk.ProjectionX().GetBinLowEdge(nbins) + hMEmultk.ProjectionX().GetBinWidth(nbins))
        for binidx in range(0, hMEmultk.ProjectionY().GetNbinsX(), binwidth):
            hSEbinmult = hSEmultk.ProjectionX(f'{comb}SEdistr', binidx, binidx + binwidth)
            hMEbinmult = hMEmultk.ProjectionX(f'{comb}MEdistr', binidx, binidx + binwidth)
            hMEreweightk.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
        hMErew[comb][region] = hMEreweightk
   
# Sum pair and antipair
for comb in combs:
    hSE[firstpair] = {}
    hSE[secondpair] = {}
    hME[firstpair] = {}
    hME[secondpair] = {}
    hMErew[firstpair] = {}
    hMErew[secondpair] = {}

    for region in regions:
        hSE[firstpair][region] = hSE[cfg['pairs'][0]][region] + hSE[cfg['pairs'][3]][region]
        hME[firstpair][region] = hME[cfg['pairs'][0]][region] + hME[cfg['pairs'][3]][region]
        hMErew[firstpair][region] = hMErew[cfg['pairs'][0]][region] + hMErew[cfg['pairs'][3]][region]
        
        hSE[secondpair][region] = hSE[cfg['pairs'][1]][region] + hSE[cfg['pairs'][2]][region]
        hME[secondpair][region] = hME[cfg['pairs'][1]][region] + hME[cfg['pairs'][2]][region]
        hMErew[secondpair][region] = hMErew[cfg['pairs'][1]][region] + hMErew[cfg['pairs'][2]][region]

# Compute the CF and write to file
for comb in list(combs.keys()) + [firstpair, secondpair]:
    for region in regions:
        rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
        print(rebin)
        hSE[comb][region].Rebin(rebin)
        hME[comb][region].Rebin(rebin)
        hMErew[comb][region].Rebin(rebin)

        if cfg['norm'] is None:  # if not specified, normalize to the yields
            norm = hMErew[comb][region].GetEntries() / hSE[comb][region].GetEntries()
        else:
            log.critical('Normalization method not implemented')

        hSE[comb][region].Sumw2()
        hCF = norm * hSE[comb][region] / hMErew[comb][region]

        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSE[comb][region].SetName('hCF')
        hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hME[comb][region].SetName('hCF')
        hME[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        
        hMErew[comb][region].SetName('hCF')
        hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hCF.SetName('hCF')
        hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCF.Write()

oFile.Close()
print(f'output saved in {oFileName}')
