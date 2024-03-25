'''
Script to compute the Lambda parameters of a pair.
The output file is LambdaParams_suffix.root

Usage:
python3 ComputeLambdaPars.py cfg.yml

'''
import os
import argparse
import yaml
import numpy as np

from ROOT import TFile, TH2D

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

# Load file containing the purities of the particles
inFilePurities = TFile(cfg['purityfilename'])

# Define the output file
oFileBaseName = 'LambdaParams'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

for pair in cfg['pairs']:

    pairCode = str(pair['code'])
    labels = ['prim', 'res', 'trueprim', 'sec']
    nBins = len(labels) 
    hPairLambdaPars = TH2D('hLambdaPars_' + pairCode, 'hLambdaPars_' + pairCode, nBins, 0, 4, nBins, 0, 4)
    hPairLambdaPars.GetXaxis().SetTitle("Particle_" + pairCode[0])
    hPairLambdaPars.GetYaxis().SetTitle("Particle_" + pairCode[1])
    for iBin, iLabel in zip(range(nBins), labels):
        hPairLambdaPars.GetXaxis().SetBinLabel(iBin+1, iLabel)
        hPairLambdaPars.GetYaxis().SetBinLabel(iBin+1, iLabel)
    
    hPurity1 = Load(inFilePurities, pair['pairpurity1'])
    purity1 = hPurity1.GetBinContent(1)
    purity1Err = hPurity1.GetBinError(1)
    hPurity2 = Load(inFilePurities, pair['pairpurity2'])
    purity2 = hPurity2.GetBinContent(1)
    purity2Err = hPurity2.GetBinError(1)
    purityFactor = purity1 * purity2
    
    primPartFracs1 = []
    primPartFracs1.append(pair['primfracs'][0])
    primPartFracs1.append(pair['resfracs'][0])
    primPartFracs1.append(pair['primfracs'][0] - pair['resfracs'][0])
    primPartFracs1.append(1-pair['primfracs'][0])

    primPartFracs2 = []
    primPartFracs2.append(pair['primfracs'][1])
    primPartFracs2.append(pair['resfracs'][1])
    primPartFracs2.append(pair['primfracs'][1] - pair['resfracs'][1])
    primPartFracs2.append(1-pair['primfracs'][1])

    
    LambdaPars = []
    LambdaParsErrs = []

    for iFrac1 in range(len(primPartFracs1)):
        for iFrac2 in range(len(primPartFracs2)):
            lambdaPar = purityFactor*primPartFracs1[iFrac1]*primPartFracs2[iFrac2]
            lambdaParErr = primPartFracs1[iFrac1]*primPartFracs2[iFrac2]*np.sqrt((purity1**2)*(purity2Err**2) + (purity2**2)*(purity1Err**2))
            hPairLambdaPars.SetBinContent(iFrac1+1, iFrac2+1, lambdaPar)
            hPairLambdaPars.SetBinError(iFrac1+1, iFrac2+1, lambdaParErr)
            LambdaPars.append(lambdaPar)
            LambdaParsErrs.append(lambdaParErr)

    hPairLambdaPars.SetStats(0)
    hPairLambdaPars.Write()

oFile.Close()
print(f'output saved in {oFileName}')
