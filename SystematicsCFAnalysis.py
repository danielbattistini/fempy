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
import numpy as np
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

# Load the CF selected for the analysis
inFileCFSelected = TFile(cfg['infileCFselected'])

# TH2D to store deviations from selected CF
systVars = list(chain.from_iterable(cfg['suffixessyst']))
nSystVars =  max(systVars) - min(systVars) + 1
nBinsKStar = len(range(cfg['systevalmaxbin']))
uppEdgeKStar = cfg['systevalmaxbin']*cfg['binwidth']

# check if the CFs for systematic variations have already been calculated, 
# else compute them
oFileCFBaseName = 'RawCF' + f'_{cfg["suffix"]}'
oFileCFName = os.path.join(f"{cfg['odir']}/systCFs", oFileCFBaseName + '_SystVarX' + '.root')
for iSystVarsAN, iANsuffix in zip(cfg['suffixessyst'], cfg['ANsuffixes']):
    for iSystVar in iSystVarsAN:
        if os.path.exists(os.path.expanduser(oFileCFName.replace('X', f"{iSystVar}"))):
            print(f'Systematics CFs for variation {iSystVar} already computed!')
        else:
            print(f'Systematics CFs for variation {iSystVar} not yet computed!')
            systDirectory = f"{cfg['odir']}/systCFs"
            try:
                os.makedirs(systDirectory, exist_ok=True)
                print(f"Directory '{systDirectory}' created successfully")
            except Exception as e:
                print(f"Error creating directory: {e}")

            os.system(f"python3 ComputeRawCF.py {args.cfg} --systvar {iSystVar} --ANsuffix {iANsuffix} ")
            print('Computed!')

inFileCFSystVar = TFile(oFileCFName.replace('X', f"{cfg['suffixessyst'][0][0]}"))
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

# Study the correlation functions for all systematic variations
regions = ["sgn", "sgn/mT0", "sgn/mT1", "sgn/mT2", "sgn/mT3", "sgn/Common", "sgn/NonCommon"]
for comb in combs:

    for region in regions:

        print(f'Checking region: {comb}/{region}/hCFrew')
        try: 
            Load(inFileCFSelected, f'{comb}/{region}/hCFrew')
            hasRegion = True
        except NameError: # Ancestors are not available, just move on
            log.warning(f"{region} not found, moving on")
            hasRegion = False

        if(hasRegion):
            print(f'Picking pair {comb}')
            oFile.mkdir(f'{comb}')
            oFile.cd(f'{comb}')

            hSystUnc = TH1D("hSystUnc", "hSystUnc", nBinsKStar, 0, uppEdgeKStar)
            hRelSystUnc = TH1D("hRelSystUnc", "hRelSystUnc", nBinsKStar, 0, uppEdgeKStar)
            hRelStatUnc = TH1D("hRelStatUnc", "hRelStatUnc", nBinsKStar, 0, uppEdgeKStar)
            hRatioStatSystUnc = TH1D("hRatioStatSystUnc", "hRatioStatSystUnc", nBinsKStar, 0, uppEdgeKStar)
            hCFStatUnc = TH1D("hCFStatUnc", "hCFStatUnc", nBinsKStar, 0, uppEdgeKStar)
            hCFSystUnc = TH1D("hCFSystUnc", "hCFSystUnc", nBinsKStar, 0, uppEdgeKStar)
            hCFStatSystUnc = TH1D("hCFStatSystUnc", "hCFStatSystUnc", nBinsKStar, 0, uppEdgeKStar)
            hResiduals = TH2D("hResiduals", "Syst vars residuals", nBinsKStar, 0, uppEdgeKStar, 
                              nSystVars, min(systVars), max(systVars))
            hResiduals.SetStats(0)

            hCFSelected = Load(inFileCFSelected, f'{comb}/{region}/hCFrew')

            if cfg.get('gaussianfits'):
                # list of histograms each containing the CF entries for
                # all variations for a specific bin 
                hCFBinEntries = []
                for iBin in range(cfg['systevalmaxbin']):
                    hCFBinEntries.append(TH1D(f"h{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV", 
                                              f"h{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV", 
                                              cfg['systhistobins'], hCFSelected.GetBinContent(iBin+1)-cfg['systhistointerval'], 
                                              hCFSelected.GetBinContent(iBin+1)+cfg['systhistointerval']))

                for iSystVar in systVars:
                    iSystVarCFFile = TFile(oFileCFName.replace('X', f"{iSystVar}"))
                    print(f'Picking syst variation number {iSystVar}')
                    hSystVar = Load(iSystVarCFFile, f'{comb}/{region}/hCFrew')
                    print(f'Picked syst variation number {iSystVar}')
                    hResiduals.GetYaxis().SetBinLabel(iSystVar-min(systVars)+1, str(iSystVar))

                    for iBin in range(nBinsKStar):
                        # the bin content of the correlation function can  
                        # be zero in the lowest mT bin at very low kstar
                        hCFBinEntries[iBin].Fill(hSystVar.GetBinContent(iBin+1))
                        hResiduals.Fill(cfg['binwidth']*(iBin+1)+cfg['binwidth']/2, iSystVar,
                                        0 if hCFSelected.GetBinContent(iBin+1) == 0 else 
                                        (hCFSelected.GetBinContent(iBin+1)-hSystVar.GetBinContent(iBin+1)) / 
                                         hCFSelected.GetBinContent(iBin+1))

                hCFYields = TH1D("hCFYields", "hCFYields", nBinsKStar, 0, uppEdgeKStar)
                hCFMeans = TH1D("hCFmeans", "hCFmeans", nBinsKStar, 0, uppEdgeKStar)

                print('Evaluating the single bins ...')
                oFile.mkdir(f'{comb}/{region}/binsfits')
                for iBin, ihCFbin in enumerate(hCFBinEntries):
                    oFile.cd(f'{comb}/{region}/binsfits')
                    iBinSelected = hCFSelected.GetBinContent(iBin+1)
                    gaus = TF1(f"gaus_{iBin}", "gaus", ihCFbin.GetBinLowEdge(1), 
                               ihCFbin.GetBinLowEdge(ihCFbin.GetNbinsX()) + ihCFbin.GetBinWidth(1))
                    gaus.SetParameter(1, iBinSelected)
                    ihCFbin.Fit(gaus, "SMRL+q", "")
                    ihCFbin.Write(f"hCFbin{round(cfg['binwidth']*iBin+(cfg['binwidth']/2))}MeV")
                    hCFYields.SetBinContent(iBin+1, gaus.GetParameter(0))
                    hCFMeans.SetBinContent(iBin+1, gaus.GetParameter(1))
                    hSystUnc.SetBinContent(iBin+1, gaus.GetParameter(2))
                    
                    # the bin content of the correlation function can  
                    # be zero in the lowest mT bin at very low kstar
                    hRelSystUnc.SetBinContent(iBin+1, 0 if iBinSelected == 0 else hSystUnc.GetBinContent(iBin+1)/iBinSelected)
                    hRelStatUnc.SetBinContent(iBin+1, 0 if iBinSelected == 0 else hCFSelected.GetBinError(iBin+1)/iBinSelected)
                    hRatioStatSystUnc.SetBinContent(iBin+1, hCFSelected.GetBinError(iBin+1)/hSystUnc.GetBinContent(iBin+1))
                    hCFStatUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFStatUnc.SetBinError(iBin+1, hCFSelected.GetBinError(iBin+1))
                    hCFSystUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFSystUnc.SetBinError(iBin+1, hSystUnc.GetBinContent(iBin+1)**2)
                    hCFStatSystUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFStatSystUnc.SetBinError(iBin+1, math.sqrt(hCFSelected.GetBinError(iBin+1)**2 + 
                                                                   hSystUnc.GetBinContent(iBin+1)**2 ) )

                hCFYields.Write()
                hCFMeans.Write()

            else: 
                binsStdDevs = []
                variedCFsEntries = []
        
                for iSystVar in systVars:
                    iSystVarBinEntries = []
                    iSystVarCFFile = TFile(oFileCFName.replace('X', f"{iSystVar}"))
                    print(f'Picking syst variation number {iSystVar}')
                    hSystVar = Load(iSystVarCFFile, f'{comb}/{region}/hCFrew')
                    print(f'Picked syst variation number {iSystVar}')
                    hResiduals.GetYaxis().SetBinLabel(iSystVar-min(systVars)+1, str(iSystVar))
                
                    for iBin in range(nBinsKStar):
                        iSystVarBinEntries.append(hSystVar.GetBinContent(iBin+1))
                        hResiduals.Fill(cfg['binwidth']*(iBin+1)+cfg['binwidth']/2, iSystVar,
                                        0 if hCFSelected.GetBinContent(iBin+1) == 0 else 
                                        (hCFSelected.GetBinContent(iBin+1)-hSystVar.GetBinContent(iBin+1)) / 
                                         hCFSelected.GetBinContent(iBin+1))

                    variedCFsEntries.append(iSystVarBinEntries)

                # print("Varied CFs entries")
                # print(variedCFsEntries)
                binsVarEntries = [[variedCFsEntries[j][i] for j in range(len(variedCFsEntries))] for i in range(len(variedCFsEntries[0]))]
                # print(binsVarEntries)
                for iBin, binVar in enumerate(binsVarEntries):
                    iBinSelected = hCFSelected.GetBinContent(iBin+1)
                    # print(f'binVar{iBin}: ')
                    # print(binVar)
                    binsStdDevs.append(np.std(binVar))
                    hSystUnc.SetBinContent(iBin+1, binsStdDevs[-1])

                    # the bin content of the correlation function can  
                    # be zero in the lowest mT bin at very low kstar
                    hRelSystUnc.SetBinContent(iBin+1, 0 if iBinSelected == 0 else hSystUnc.GetBinContent(iBin+1)/iBinSelected)
                    hRelStatUnc.SetBinContent(iBin+1, 0 if iBinSelected == 0 else hCFSelected.GetBinError(iBin+1)/iBinSelected)
                    hRatioStatSystUnc.SetBinContent(iBin+1, 0 if iBinSelected == 0 else hCFSelected.GetBinError(iBin+1)/hSystUnc.GetBinContent(iBin+1))
                    hCFStatUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFStatUnc.SetBinError(iBin+1, hCFSelected.GetBinError(iBin+1))
                    hCFSystUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFSystUnc.SetBinError(iBin+1, hSystUnc.GetBinContent(iBin+1))
                    hCFStatSystUnc.SetBinContent(iBin+1, iBinSelected)
                    hCFStatSystUnc.SetBinError(iBin+1, math.sqrt(hCFSelected.GetBinError(iBin+1)**2 + 
                                                                     hSystUnc.GetBinContent(iBin+1)**2 ) )

            oFile.mkdir(f'{comb}/{region}/unc')
            oFile.cd(f'{comb}/{region}/unc')
            print('Writing histos with error informations ...')
            hSystUnc.Write()
            hRelSystUnc.Write()
            hRelStatUnc.Write()
            hRatioStatSystUnc.Write()
            hCFStatUnc.Write()
            hCFSystUnc.Write()
            hCFStatSystUnc.Write()
            hResiduals.Write()

# oFile.mkdir('compare_residuals')
# oFile.cd('compare_residuals')
# for iComb in range(0, len(combs), 2):
#     hResiduals = TH2D(f"hRatioResiduals_{combs[iComb]}_{combs[iComb+1]}", 
#                       "Syst vars residuals", nBinsKStar, 0, uppEdgeKStar, 
#                       nSystVars, min(systVars), max(systVars))
#     hResiduals.SetStats(0)
#     hResidualsPairOne = Load(oFile, f"{combs[iComb]}/unc/hResiduals")
#     hResidualsPairTwo = Load(oFile, f"{combs[iComb+1]}/unc/hResiduals")
#     for iBinX in range(nBinsKStar):
#         for iBinY in range(nSystVars):
#             resCombOne = hResidualsPairOne.GetBinContent(iBinX+1, iBinY+1)
#             resCombTwo = hResidualsPairTwo.GetBinContent(iBinX+1, iBinY+1)
#             if resCombOne!=0 and resCombOne!=0: 
#                 hResiduals.SetBinContent(iBinX+1, iBinY+1, resCombOne/resCombTwo)
#             else:
#                 hResiduals.SetBinContent(iBinX+1, iBinY+1, 0)
#     hResiduals.Write()

oFile.Close()
print(f'output saved in {oFileName}')
