'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import argparse
import yaml

from ROOT import TFile, TH2F, TH1D

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

regions = ['CFsystvar']

# Define the output file
oFileBaseName = 'SystCFs'
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

# list of histograms each containing the CF entries for
# all variations for a specific bin 
hCFBinEntries = []
for iBin in range(cfg['systevalnbins']):
    hCFBinEntries.append(TH1D(f"h{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV", 
                              f"h{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV", 
                              cfg['systhistobins'], cfg['systhistolowedge'], 
                              cfg['systhistouppedge']))

# Compute the correlation functions for all systematic variations
hSE = {}
hME = {}
hMErew = {}
for iSystFile, systFileName in enumerate(cfg['infilessyst']):
    print(f'Picking file {iSystFile}')
    systFile = TFile(systFileName)
    for iSystVar in cfg['suffixessyst'][iSystFile]:
        print(f'runsuffix {iSystVar}')
        for comb in combs:
            iPart1 = int(comb[1])
            iPart2 = int(comb[2])
            oFile.cd(comb)

            hSE[comb] = {}
            hME[comb] = {}
            hMErew[comb] = {}
            for region in regions:
                if f'HMResults{iSystVar}' in GetKeyNames(systFile): # Make correlation functions from FemtoDream
                    fdcomb = f'Particle{iPart1}_Particle{iPart2}'
                    # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
                    hSE[comb][region] = TH1D()
                    Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/SEDist_{fdcomb}').Copy(hSE[comb][region])

                    # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
                    hME[comb][region] = TH1D()
                    Load(systFile, f'HMResults{iSystVar}/HMResults{iSystVar}/{fdcomb}/MEDist_{fdcomb}').Copy(hME[comb][region])

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
        for comb in combs + ['p02_13', 'p03_12']:
            for region in regions:
                rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
                hSE[comb][region].Rebin(rebin)
                hME[comb][region].Rebin(rebin)
                hMErew[comb][region].Rebin(rebin)

                if cfg['norm'] is None:  # if not specified, normalize to the yields
                    norm = hME[comb][region].GetEntries() / hSE[comb][region].GetEntries()
                    normrew = hMErew[comb][region].GetEntries() / hSE[comb][region].GetEntries()
                else:
                    firstBin = hSE[comb][region].FindBin(cfg['norm'][0]*1.0001)
                    lastBin = hSE[comb][region].FindBin(cfg['norm'][1]*0.9999)
                    norm = hME[comb][region].Integral(firstBin, lastBin) / hSE[comb][region].Integral(firstBin, lastBin)
                    normrew = hMErew[comb][region].Integral(firstBin, lastBin) / hSE[comb][region].Integral(firstBin, lastBin)

                hSE[comb][region].Sumw2()
                hME[comb][region].Sumw2() # Just to trigger the same draw option as for hSE

                hCF = norm * hSE[comb][region] / hME[comb][region]
                hCFrew = normrew * hSE[comb][region] / hMErew[comb][region]

                oFile.mkdir(f'{comb}/{region}/{iSystVar}')
                oFile.cd(f'{comb}/{region}/{iSystVar}')

                hSE[comb][region].SetName(f'hSE_{iSystVar}')
                hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
                hSE[comb][region].Write()

                hME[comb][region].SetName(f'hME_{iSystVar}')
                hME[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
                hME[comb][region].Write()

                hMErew[comb][region].SetName(f'hMErew_{iSystVar}')
                hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

                hCF.SetName(f'hCF_{iSystVar}')
                hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
                hCF.Write()

                hCFrew.SetName(f'hCFrew_{iSystVar}')
                hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
                hCFrew.Write()

                for iBin in range(cfg['systevalnbins']):
                    hCFBinEntries[iBin].Fill(hCFrew.GetBinContent(iBin+1))

                print(f'Finished systematic variation n° {iSystVar}')

#hCFselected = TH1D() Load(TFile(cfg['infileCF']), cfg['infileCFpath'])
#
#oFile.cd()
#oFile.mkdir('binsfits')
#hCFYields = TH1D("hCFYields", "hCFYields", cfg['systevalnbins'], 
#                 cfg['systevalnbins'][0]*cfg['binwidth'], 
#                 cfg['systevalnbins'][1]*cfg['binwidth'])
#hCFmeans = TH1D("hCFmeans", "hCFmeans", cfg['systevalnbins'], 
#                 cfg['systevalnbins'][0]*cfg['binwidth'], 
#                 cfg['systevalnbins'][1]*cfg['binwidth'])
#hSystUnc = TH1D("hSystUnc", "hSystUnc", cfg['systevalnbins'], 
#                 cfg['systevalnbins'][0]*cfg['binwidth'], 
#                 cfg['systevalnbins'][1]*cfg['binwidth'])
#hRelSystUnc = TH1D("hRelSystUnc", "hRelSystUnc", cfg['systevalnbins'], 
#                 cfg['systevalnbins'][0]*cfg['binwidth'], 
#                 cfg['systevalnbins'][1]*cfg['binwidth'])
#hRelStatUnc = TH1D("hRelStatUnc", "hRelStatUnc", cfg['systevalnbins'], 
#                   cfg['systevalnbins'][0]*cfg['binwidth'], 
#                   cfg['systevalnbins'][1]*cfg['binwidth'])
#hRatioStatSystUnc = TH1D("hRatioStatSystUnc", "hRatioStatSystUnc", cfg['systevalnbins'], 
#                         cfg['systevalnbins'][0]*cfg['binwidth'], 
#                         cfg['systevalnbins'][1]*cfg['binwidth'])
#hCFWithStatSystUnc = TH1D("hCFWithStatSystUnc", "hCFWithStatSystUnc", cfg['systevalnbins'], 
#                          cfg['systevalnbins'][0]*cfg['binwidth'], 
#                          cfg['systevalnbins'][1]*cfg['binwidth'])
#for iBin, ihCFbin in enumerate(hCFBinEntries):
#    oFile.cd('binsfits')
#    ihCFbin.Fit("gaus")
#    ihCFbin.Write(f"hCFbin{cfg['binwidth']*iBin+(cfg['binwidth']/2)}MeV")
#    gausFit = ihCFbin.GetListOfFunctions().FindObject("gaus")
#    hCFYields.SetBinContent(iBin+1, gausFit.GetParameter(0))
#    hCFMeans.SetBinContent(iBin+1, gausFit.GetParameter(1))
#    hSystUnc.SetBinContent(iBin+1, gausFit.GetParameter(2))
#    hRelSystUnc.SetBinContent(iBin+1, hSystUnc.GetBinContent(iBin+1)/hCFselected.GetBinContent(iBin+1))
#    hRelStatUnc.SetBinContent(iBin+1, hCFselected.GetBinError(iBin+1)/hCFselected.GetBinContent(iBin+1))
#    hRelStatUnc.SetBinContent(iBin+1, hCFselected.GetBinError(iBin+1)/hSystUnc.GetBinContent(iBin+1))
#    hCFWithStatSystUnc.SetBinContent(iBin+1, hCFselected.GetBinContent(iBin+1))
#    hCFWithStatSystUnc.SetBinError(iBin+1, np.sqrt(hCFselected.GetBinError(iBin+1)**2 + 
#                                                   hSystUnc.GetBinContent(iBin+1)**2 ) )
#    
#oFile.mkdir('systunc')
#oFile.cd('systunc')    
#hCFYields.Write()
#hCFMeans.Write()
#hSystUnc.Write()
#hRelSystUnc.Write()
#hRelStatUnc.Write()
#hRelStatUnc.Write()
#hCFWithStatSystUnc.Write()
#
oFile.Close()
print(f'output saved in {oFileName}')
