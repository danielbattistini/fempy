'''
Script to compute the raw correlation function.
The output file is RawSystCF_suffix.root

Usage:
python3 ComputeSystematicsCF.py cfg.yml

'''
import os
import argparse
import yaml
from itertools import chain

from ROOT import TFile, TH2F, TH1D, TH2D

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
oFileBaseName = 'RawSystCF'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)
combs = ['p02', 'p03', 'p12', 'p13']
for comb in combs:
    oFile.mkdir(comb)

totalSystVars = len(list(chain.from_iterable(cfg['suffixessyst'])))
hFemtoPairs = TH2D("hFemtoPairs", "Pairs in femto region", totalSystVars, 0, totalSystVars, len(combs), 0, len(combs))
hFemtoPairs.SetStats(0)
hFemtoPairs.GetYaxis().SetLabelSize(50)
hFemtoPairs.GetXaxis().SetLabelSize(30)
hFemtoPairs.GetXaxis().CenterLabels()

# Compute the correlation functions for all systematic variations
hSE = {}
hME = {}
hMErew = {}
for iSystFile, systFileName in enumerate(cfg['infilessyst']):
    systFile = TFile(systFileName)
    for iSystVar in cfg['suffixessyst'][iSystFile]:
        for iComb, comb in enumerate(combs):
            iPart1 = int(comb[1])
            iPart2 = int(comb[2])
            oFile.cd(comb)

            if f'HMResults{iSystVar}' in GetKeyNames(systFile): # Make correlation functions from FemtoDream
                fdcomb = f'Particle{iPart1}_Particle{iPart2}'
                # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
                hSE[comb] = TH1D()
                Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/SEDist_{fdcomb}').Copy(hSE[comb])


                hFemtoPairs.Fill(iSystVar-0.5, len(combs)-iComb-0.5, 
                                 hSE[comb].Integral(hSE[comb].FindBin(0.0001), hSE[comb].FindBin(0.2*0.9999)))
                hFemtoPairs.GetYaxis().SetBinLabel(len(combs)-iComb, comb)
                hFemtoPairs.GetXaxis().SetBinLabel(iSystVar, str(iSystVar))

                # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
                hME[comb] = TH1D()
                Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/MEDist_{fdcomb}').Copy(hME[comb])

                # No need to cast the these because the projection of TH2D and TH2F is always TH1D
                hSEmultk = Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/SEMultDist_{fdcomb}')
                hMEmultk = Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/MEMultDist_{fdcomb}')

            nbins = hMEmultk.ProjectionX().GetNbinsX()
            hMEreweightk = TH1D("MErewdistr", "MErewdistr", nbins, hMEmultk.GetXaxis().GetXmin(), hMEmultk.GetXaxis().GetXmax())
            for iBin in range(hMEmultk.ProjectionY().GetNbinsX() + 2): # Loop over underflow, all bins, and overflow
                hSEbinmult = hSEmultk.ProjectionX(f'{comb}SEdistr', iBin, iBin)
                hMEbinmult = hMEmultk.ProjectionX(f'{comb}MEdistr', iBin, iBin)
                if(hMEbinmult.Integral() > 0):
                    hMEreweightk.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
            hMErew[comb] = hMEreweightk

        # Sum pair and antipair
        for comb in combs:
            hSE['p02_13'] = hSE['p02'] + hSE['p13']
            hME['p02_13'] = hME['p02'] + hME['p13']
            hMErew['p02_13'] = hMErew['p02'] + hMErew['p13']

            hSE['p03_12'] = hSE['p03'] + hSE['p12']
            hME['p03_12'] = hME['p03'] + hME['p12']
            hMErew['p03_12'] = hMErew['p03'] + hMErew['p12']

        # Compute the CF and write to file
        for comb in combs + ['p02_13', 'p03_12']:
            rebin = round(float(cfg['binwidth']) / (hSE[comb].GetBinWidth(1) * 1000))
            hSE[comb].Rebin(rebin)
            hME[comb].Rebin(rebin)
            hMErew[comb].Rebin(rebin)

            if cfg['norm'] is None:  # if not specified, normalize to the yields
                norm = hME[comb].GetEntries() / hSE[comb].GetEntries()
                normrew = hMErew[comb].GetEntries() / hSE[comb].GetEntries()
            else:
                firstBin = hSE[comb].FindBin(cfg['norm'][0]*1.0001)
                lastBin = hSE[comb].FindBin(cfg['norm'][1]*0.9999)
                norm = hME[comb].Integral(firstBin, lastBin) / hSE[comb].Integral(firstBin, lastBin)
                normrew = hMErew[comb].Integral(firstBin, lastBin) / hSE[comb].Integral(firstBin, lastBin)

            hSE[comb].Sumw2()
            hME[comb].Sumw2() # Just to trigger the same draw option as for hSE

            hCF = norm * hSE[comb] / hME[comb]
            hCFrew = normrew * hSE[comb] / hMErew[comb]

            oFile.mkdir(f'{comb}/var{iSystVar}')
            oFile.cd(f'{comb}/var{iSystVar}')

            hSE[comb].SetName(f'hSE_{iSystVar}')
            hSE[comb].SetTitle(';#it{k}* (GeV/#it{c});Counts')
            hSE[comb].Write()

            hME[comb].SetName(f'hME_{iSystVar}')
            hME[comb].SetTitle(';#it{k}* (GeV/#it{c});Counts')
            hME[comb].Write()

            hMErew[comb].SetName(f'hMErew_{iSystVar}')
            hMErew[comb].SetTitle(';#it{k}* (GeV/#it{c});Counts')

            hCF.SetName(f'hCF_{iSystVar}')
            hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
            hCF.Write()

            hCFrew.SetName(f'hCFrew_{iSystVar}')
            hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
            hCFrew.Write()

        print(f'Finished systematic variation {iSystVar}')

oFile.cd()
hFemtoPairs.Write()
oFile.Close()
print(f'output saved in {oFileName}')
