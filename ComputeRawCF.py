'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import argparse
import yaml

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

regions = ['sgn']
combs = ['p02', 'p03', 'p12', 'p13']

# Load input file with same- and mixed-event distributions
inFile = TFile(cfg['infile'])
runSuffix = cfg['runsuffix']

try: # Check if there are ancestor histograms
    # Only load SE, the ME is the same for all
    Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/Particle0_Particle2/SEDistCommon_Particle0_Particle2')
    hasAncestors = True
    ancestors = ['Common', 'NonCommon']
except NameError: # Ancestors are not available, just move on
    log.warning("Ancestors not found, moving on")
    hasAncestors = False
    ancestors = []

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
hWeightsRew = {}
hSEmultk = {}
hMEmultk = {}

# 5 bins: sgn, sgn_Common, sgn_NonCommon, sbl, sbr. If some are not present they are left empty
hFemtoPairs = TH2D("hFemtoPairs", "Pairs in femto region", 5, 0, 5, len(combs) + int(len(combs)/2), 0, len(combs) + int(len(combs)/2))
hFemtoPairs.GetXaxis().SetBinLabel(1, "sgn")
hFemtoPairs.GetXaxis().SetBinLabel(2, "sgn/Common")
hFemtoPairs.GetXaxis().SetBinLabel(3, "sgn/NonCommon")
hFemtoPairs.GetXaxis().SetBinLabel(4, "sbl")
hFemtoPairs.GetXaxis().SetBinLabel(5, "sbr")
hFemtoPairs.GetYaxis().SetLabelSize(50)
hFemtoPairs.GetXaxis().SetLabelSize(50)

for ncomb, comb in enumerate(combs):
    iPart1 = int(comb[1])
    iPart2 = int(comb[2])
    oFile.mkdir(comb)
    oFile.cd(comb)

    hSE[comb] = {}
    hME[comb] = {}
    hMErew[comb] = {}
    hWeightsRew[comb] = {}
    hSEmultk[comb] = {}
    hMEmultk[comb] = {}

    # Load the same- and mixed-event distributions for each region (signal/sidebands)
    for iRegion, region in enumerate(regions):
        runSuffix = cfg['runsuffix']
        if f'HMResults{runSuffix}' in GetKeyNames(inFile): # Make correlation functions from FemtoDream
            fdcomb = f'Particle{iPart1}_Particle{iPart2}'
            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hME[comb][region] = TH1D()
            Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/MEDist_{fdcomb}').Copy(hME[comb][region])
            hMEmultk[comb] = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/MEMultDist_{fdcomb}')
            hSE[comb][region] = TH1D()
            Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEDist_{fdcomb}').Copy(hSE[comb][region])

            if 'multrewMCwithdata' in cfg:
                log.info('Reweighting MC with data!')
                inFileData = TFile(cfg['infilerew'])
                runSuffixData = cfg['runsuffixrew']
                # Only load the SE vs mult, which is needed for the weights
                hSEmultk[comb] = Load(inFileData, f'HMResults{runSuffixData}/HMResults{runSuffixData}/{fdcomb}/SEMultDist_{fdcomb}')
                hSEmultk[comb].SetDirectory(0)
                inFileData.Close()
            else:
                hSEmultk[comb] = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEMultDist_{fdcomb}')
                hSEmultk[comb].SetDirectory(0)
            # Do ancestors CF if they are available
            for ancestor in ancestors:
                # Only load SE, the ME is the same as before
                hSEAnc = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEDist{ancestor}_{fdcomb}')
                hSEAncMultk = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEMultDist{ancestor}_{fdcomb}')
                hSE[comb][f'{region}/{ancestor}'] = TH1D()
                hSEmultk[comb][f'{region}/{ancestor}'] = TH1D()
                hSEAnc.Copy(hSE[comb][f'{region}/{ancestor}'])
                hSEAncMultk.Copy(hSEmultk[comb][f'{region}/{ancestor}'])

        elif comb in GetKeyNames(inFile): # Make correlation functions from ALICE3 simulations
            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hSE[comb][region] = TH1D()
            Load(inFile, f'{comb}/hSE').Copy(hSE[comb][region])

            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hME[comb][region] = TH1D()
            Load(inFile, f'{comb}/hME').Copy(hME[comb][region])

            # Mult reweighting not implemented. Keep dummy histogram for now
            nbins = hSE[comb][region].GetNbinsX()
            xMin = hSE[comb][region].GetXaxis().GetXmin()
            xMax = hSE[comb][region].GetXaxis().GetXmax()
            hSEmultk[comb] = TH2F('hSEMult', '', nbins, xMin, xMax, 200, 0, 200)
            hMEmultk[comb] = TH2F('hMEMult', '', nbins, xMin, xMax, 200, 0, 200)

        nbins = hMEmultk[comb].ProjectionX().GetNbinsX()
        hWeights = TH1D(f"Weights{region}_{comb}", f"Weights{region}_{comb}", hMEmultk[comb].ProjectionY().GetNbinsX() + 2, 
                        hMEmultk[comb].GetYaxis().GetXmin(), hMEmultk[comb].GetYaxis().GetXmax())
        hMEreweightk = TH1D(f"MErewdistr{region}_{comb}", f"MErewdistr{region}_{comb}", nbins, hMEmultk[comb].GetXaxis().GetXmin(), hMEmultk[comb].GetXaxis().GetXmax())

        startBinMultRew, endBinMultRew = cfg.get('binmultrew', [0, hMEmultk[comb].ProjectionY().GetNbinsX() + 2]) 
        for iBin in range(startBinMultRew, endBinMultRew): # Loop over underflow, all bins, and overflow
            hSEbinmult = hSEmultk[comb].ProjectionX(f'{comb}SEdistr', iBin, iBin)
            hMEbinmult = hMEmultk[comb].ProjectionX(f'{comb}MEdistr', iBin, iBin)

            if hMEbinmult.Integral() > 0:
                hMEreweightk.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
                hWeights.SetBinContent(iBin, hSEbinmult.Integral()/hMEbinmult.Integral())

                # Save the CFs for each multiplicity bin
                oFile.mkdir(f'{comb}/multbins/{iBin}')
                oFile.cd(f'{comb}/multbins/{iBin}')
                firstBin = hSEbinmult.FindBin(cfg['norm'][0] * 1.0001)
                lastBin = hSEbinmult.FindBin(cfg['norm'][1] * 0.9999)
                norm = hMEbinmult.Integral(firstBin, lastBin) / hSEbinmult.Integral(firstBin, lastBin)
                hCFbinmult = norm * hSEbinmult / hMEbinmult
                hCFbinmult.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
                hCFbinmult.Write(f'hCF_multbin{iBin}')
                hSEbinmult.Write(f'hSE_multbin{iBin}')
                hMEbinmult.Write(f'hME_multbin{iBin}')
                oFile.cd(comb)

        hMErew[comb][region] = hMEreweightk
        hWeightsRew[comb][region] = hWeights

# Sum pair and antipair
for comb in combs:
    hSE['p02_13'] = {}
    hSE['p03_12'] = {}
    hME['p02_13'] = {}
    hME['p03_12'] = {}
    hMErew['p02_13'] = {}
    hMErew['p03_12'] = {}
    hWeightsRew['p02_13'] = {}
    hWeightsRew['p03_12'] = {}

    for region in regions:
        hSE['p02_13'][region] = hSE['p02'][region] + hSE['p13'][region]
        hME['p02_13'][region] = hME['p02'][region] + hME['p13'][region]
        hMErew['p02_13'][region] = hMErew['p02'][region] + hMErew['p13'][region]
        hWeightsRew['p02_13'][region] = hWeightsRew['p02'][region] + hWeightsRew['p13'][region]
        
        hSE['p03_12'][region] = hSE['p03'][region] + hSE['p12'][region]
        hME['p03_12'][region] = hME['p03'][region] + hME['p12'][region]
        hMErew['p03_12'][region] = hMErew['p03'][region] + hMErew['p12'][region]
        hWeightsRew['p03_12'][region] = hWeightsRew['p03'][region] + hWeightsRew['p12'][region]

# Compute the CF and write to file
for iComb, comb in enumerate(combs + ['p02_13', 'p03_12']):
    for region in hSE[comb]:
        rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
        hSE[comb][region].Rebin(rebin)
        regionME = region.split('/')[0] # remove '/Common' and '/NonCommon' since the ME is the same as the inclusive ancestor
        hME[comb][regionME].Rebin(rebin)
        hMErew[comb][region].Rebin(rebin)

        # Count number of pairs in the femto region
        regionBin = hFemtoPairs.GetXaxis().FindBin(region)
        firstBin = hSE[comb][region].FindBin(0.0001)
        lastBin = hSE[comb][region].FindBin(0.2*0.9999)
        hFemtoPairs.GetYaxis().SetBinLabel(iComb + 1, comb)
        hFemtoPairs.SetBinContent(regionBin + 1, iComb + 1, hSE[comb][region].Integral(firstBin, lastBin))

        if cfg['norm'] is None:  # if not specified, normalize to the yields
            firstBinNorm = 1
            lastBinNorm = hMErew[comb][region].GetNbinsX()
        elif isinstance(cfg['norm'], list):
            firstBinNorm = hSE[comb][region].FindBin(cfg['norm'][0]*1.0001)
            lastBinNorm = hSE[comb][region].FindBin(cfg['norm'][1]*0.9999)
        else:
            log.critical('Normalization method not implemented')

        YieldSE = hSE[comb][region].Integral(firstBinNorm, lastBinNorm)
        norm = hME[comb][regionME].Integral(firstBinNorm, lastBinNorm) / YieldSE
        normRew = hMErew[comb][region].Integral(firstBinNorm, lastBinNorm) / YieldSE

        hSE[comb][region].Sumw2()
        hME[comb][regionME].Sumw2() # Just to trigger the same draw option as for hSE

        hCF = norm * hSE[comb][region] / hME[comb][region]
        hCFrew = normRew * hSE[comb][region] / hMErew[comb][region]

        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSE[comb][region].SetName('hSE')
        hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        hSE[comb][region].Write()

        hME[comb][regionME].SetName('hME')
        hME[comb][regionME].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        hME[comb][regionME].Write()

        hMErew[comb][region].SetName('hMErew')
        hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hWeightsRew[comb][region].SetName('hWeightsRew')
        hWeightsRew[comb][region].SetTitle(';Mult bin;Weight')
        hWeightsRew[comb][region].Write()

        hCF.SetName('hCF')
        hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCF.Write()

        hCFrew.SetName('hCFrew')
        hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCFrew.Write()

oFile.cd()
hFemtoPairs.Write()
oFile.Close()
print(f'output saved in {oFileName}')

