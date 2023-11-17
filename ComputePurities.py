'''
Script to compute the purities of correlated particles.
The output file is RawPurities_suffix.root

Usage:
python3 ComputePurities.py cfg.yml

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

# Load input file with same- and mixed-event distributions
inFileData = TFile(cfg['infiledata'])
inFileMC = TFile(cfg['infilemc'])

# Define the output file
oFileBaseName = 'RawPurities'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

runSuffix = cfg['runsuffix']

chargedpart = {
    'p0': ['Particle0',''],
    'p1': ['Particle1','Anti']
}

neutralpart = {
    'p2': ['Particle2','v0'],
    'p3': ['Particle3','Antiv0']
}

# charged particles purity
for cpart, fdcpart in chargedpart.items():
    oFile.mkdir(cpart)
    oFile.cd(cpart)
    
    hIdentPartPt = Load(inFileMC, f'HM{fdcpart[1]}TrkCutsMC{runSuffix}/HM{fdcpart[1]}TrkCutsMC{runSuffix}/IdentPartPt')
    hCorrPartPt = Load(inFileMC, f'HM{fdcpart[1]}TrkCutsMC{runSuffix}/HM{fdcpart[1]}TrkCutsMC{runSuffix}/CorrParPt')
    hCorrPartPt.Divide(hIdentPartPt)
    hIDratio = hCorrPartPt
    hIDratio.SetTitle('CorrectIdentifRatio')
    hIDratio.SetName('CorrectIdentifRatio')
    hIDratio.Write()
    for npart, fdnpart in neutralpart.items():
        oFile.mkdir(f'{cpart}/{npart}')
        oFile.cd(f'{cpart}/{npart}')
        kstarPt = Load(inFileMC, f'HMResultsQA{runSuffix}/HMResultsQA{runSuffix}/PairQA/QA_{fdcpart[0]}_{fdnpart[0]}/KstarPtSEPartOne_{fdcpart[0]}_{fdnpart[0]}')
        lowbin = kstarPt.ProjectionX().FindBin(0)
        uppbin = kstarPt.ProjectionX().FindBin(0.2)
        hpTdistr = kstarPt.ProjectionY('femto{fdcpart[0]}', lowbin, uppbin)
        hpTdistr.Scale(1/hpTdistr.Integral())
        hpTdistr.SetTitle(f'pTdistr_0_200_{fdcpart[0]}')
        hpTdistr.SetName(f'pTdistr_0_200_{fdcpart[0]}')
        hpTdistr.Write()
        weights = []
        purity = 0
        for iBin in range(hIDratio.GetNbinsX()):
            weight = hpTdistr.Interpolate(hIDratio.GetBinCenter(iBin))
            weights.append(weight)
            purity += weight * hIDratio.GetBinContent(iBin)
        purity = purity/sum(weights) 
        hPurity = TH1F('Purity', 'Purity', 1, 0, 1)
        hPurity.Fill(0, purity)
        hPurity.Write()

oFile.Close()
print(f'output saved in {oFileName}')

