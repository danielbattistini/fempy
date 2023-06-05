'''
Script to create a tree dataset from binned dataset
'''
import argparse

import numpy as np

from ROOT import TFile, TTree, gRandom

from utils.io import GetObjectFromFile
from utils.handle import GetFemtoDreamPairId


def CreateFakeDataSet(reducedDataset=False, **kwargs):
    """
    Create a fake tree dataset using a histogram as input.

    Parameters
    ----------
    inFileName : str
        path of input datafile
    inHistPath : str
        path of histogram to sample inside the input file
    oFileName : str
        path of output file
    oTreeName : str
        name of the output tree
    reducedDataset : bool
        reduce by a factor 1000 the size of datasample
    Returns
    -------
    int
        Description of return value

    """
    reducedDataset = True
    scaleFactor = 0.1

    inFile = TFile(kwargs['inFileName'])

    pairCombs = GetFemtoDreamPairId('sc')

    regions = {'sgn': '', 'sbl': 'SBLeft_', 'sbr': 'SBRight_'}

    hSE = {}
    hME = {}
    for regKey, regVal in regions.items():
        directory = f'HM_CharmFemto_{regVal}Dpion_Results0'

        # load histograms
        for pairKey, pairComb in pairCombs.items():
            print(f'{directory}/{directory}/{pairComb}/SEDist_{pairComb}')
            hSE[f'{regKey}_{pairComb}'] = GetObjectFromFile(inFile, f'{directory}/{directory}/{pairComb}/SEDist_{pairComb}')
            hME[f'{regKey}_{pairComb}'] = GetObjectFromFile(inFile, f'{directory}/{directory}/{pairComb}/MEDist_{pairComb}')
    inFile.Close()

    oFile = TFile(kwargs['oFileName'], 'recreate')
    for regKey, regVal in regions.items():
        directory = f'HM_CharmFemto_{regVal}Dpion_Results0'

        # Generate se and me distributions
        kStar = np.empty((1), dtype='d')
        for pairKey, pairComb in pairCombs.items():
            # sample SE
            nSamplSE = gRandom.Poisson(hSE[f'{regKey}_{pairComb}'].GetEntries())
            if reducedDataset:
                nSamplSE = int(nSamplSE * scaleFactor)
            print(regKey, pairKey, nSamplSE, hSE[f'{regKey}_{pairComb}'].GetEntries())
            treeSE = TTree(f'treeDpi_SE_{regKey}_{pairComb}', f'treeDpi SE {regKey} {pairComb}')
            treeSE.Branch('kStar', kStar, 'kStar/D')

            for _ in range(nSamplSE):
                kStar[0] = hSE[f'{regKey}_{pairComb}'].GetRandom() * 1000  # transorm to MeV/c
                treeSE.Fill()
            treeSE.Write()

            ## bugged
            # # sample ME
            # nSamplME = gRandom.Poisson(hME[f'{regKey}_{pairComb}'].GetEntries())
            # if reducedDataset:
            #     nSamplME = int(nSamplME * scaleFactor0)
            # print(nSamplME, hME[f'{regKey}_{pairComb}'].GetEntries())
            # treeME = TTree(f'treeDpi_ME_{regKey}_{pairComb}', f'treeDpi ME {regKey} {pairComb}')
            # treeME.Branch('kStar', kStar, 'kStar/D')

            # for _ in range(nSamplME):
            #     kStar[0] = hME[f'{regKey}_{pairComb}'].GetRandom() * 1000  # transorm to MeV/c
            #     treeME.Fill()

            # # Write to file
            # treeME.Write()

            nSamplME = gRandom.Poisson(hME[f'{regKey}_{pairComb}'].GetEntries())
            if reducedDataset:
                nSamplME = int(nSamplME * scaleFactor)
            print(regKey, pairKey, nSamplME, hME[f'{regKey}_{pairComb}'].GetEntries())
            treeME = TTree(f'treeDpi_ME_{regKey}_{pairComb}', f'treeDpi ME {regKey} {pairComb}')
            treeME.Branch('kStar', kStar, 'kStar/D')

            for _ in range(nSamplME):
                kStar[0] = hME[f'{regKey}_{pairComb}'].GetRandom() * 1000  # transorm to MeV/c
                treeME.Fill()
            treeME.Write()

    oFile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('inFileName', metavar='text', help='AnalysisResults.root file')
    parser.add_argument('oFileName', metavar='text', help='AnalysisResults.root file')
    parser.add_argument('-r', action='store_true', default=False, help='use reduced size of dataset')

    args = parser.parse_args()

    CreateFakeDataSet(inFileName=args.inFileName, oFileName=args.oFileName, reducedDataset=args.r)
