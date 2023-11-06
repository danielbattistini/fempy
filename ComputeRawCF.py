'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import argparse
import yaml

from ROOT import TFile

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
    'pp': 'Particle0_Particle2',
    'pm': 'Particle0_Particle3',
    'mp': 'Particle1_Particle2',
    'mm': 'Particle1_Particle3'
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
for comb, fdcomb in combs.items():
    oFile.mkdir(comb)
    oFile.cd(comb)
    
    hSE[comb] = {}
    hME[comb] = {}
    for region in regions:
        if(cfg['domultrew']):
            binwidth = cfg['multbinwidth']
            hSEmultk = Load(inFile, f'HMResults0/HMResults0/{fdcomb}/SEMultDist_{fdcomb}')
            hSE[comb][region] = hSEmultk.ProjectionX()
            hMEmultk = Load(inFile, f'HMResults0/HMResults0/{fdcomb}/MEMultDist_{fdcomb}')
            hMEreweightk = (hMEmultk.ProjectionX()).Clone()
            hMEreweightk.Reset('ICESM')
            for binidx in range(0, (hMEmultk.ProjectionY()).GetNbinsX(), binwidth):
                hSEbinmult = hSEmultk.ProjectionX(comb + 'SEdistr', binidx, binidx + binwidth)
                hMEbinmult = hMEmultk.ProjectionX(comb + 'MEdistr', binidx, binidx + binwidth)
                hMEbinmult.Scale(hSEbinmult.Integral()/hMEbinmult.Integral())
                hMEreweightk.Add(hMEbinmult)
            hME[comb][region] = hMEreweightk
        else:
            hSE[comb][region] = Load(inFile, f'HMResults0/HMResults0/{fdcomb}/SEDist_{fdcomb}')
            hME[comb][region] = Load(inFile, f'HMResults0/HMResults0/{fdcomb}/MEDist_{fdcomb}')

# Sum pair and antipair
for comb in combs:
    hSE['sc'] = {}
    hSE['oc'] = {}
    hME['sc'] = {}
    hME['oc'] = {}

    for region in regions:
        hSE['sc'][region] = hSE['pp'][region] + hSE['mm'][region]
        hME['sc'][region] = hME['pp'][region] + hME['pp'][region]

        hSE['oc'][region] = hSE['pm'][region] + hSE['mp'][region]
        hME['oc'][region] = hME['pm'][region] + hME['mp'][region]

# Compute the CF and write to file
for comb in list(combs.keys()) + ['sc', 'oc']:
    for region in regions:
        rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
        print(rebin)
        hSE[comb][region].Rebin(rebin)
        hME[comb][region].Rebin(rebin)

        if cfg['norm'] is None:  # if not specified, normalize to the yields
            norm = hME[comb][region].GetEntries() / hSE[comb][region].GetEntries()
        else:
            log.critical('Normalization method not implemented')

        hSE[comb][region].Sumw2()
        hCF = norm * hSE[comb][region] / hME[comb][region]

        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSE[comb][region].SetName('hCF')
        hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hME[comb][region].SetName('hCF')
        hME[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hCF.SetName('hCF')
        hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*})')
        hCF.Write()

oFile.Close()
print(f'output saved in {oFileName}')
