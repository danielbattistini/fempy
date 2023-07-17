import argparse
import os

from ROOT import TFile

import fempy
from fempy.utils.analysis import GetNormFactor

parser = argparse.ArgumentParser()
parser.add_argument('inFileName')
parser.add_argument('oFileName')
parser.add_argument('--bw', default=50, type=int)
parser.add_argument('--norm', default=(1500, 2000), type=float, nargs=2)
args = parser.parse_args()

inFile = TFile(args.inFileName)
oFile = TFile(args.oFileName, 'create')

combs = fempy.utils.io.GetCombs(inFile)
regions = fempy.utils.io.GetRegions(inFile.Get(f'{combs[0]}/SE'))

print('found the following combinations:', combs)
print('found the following regions:', regions)

bw = args.bw

for comb in combs:
    oFile.mkdir(comb)
    oFile.cd(comb)

    for region in regions:
        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSEMultVsKStar = inFile.Get(f'{comb}/SE/{region}/hCharmPtVsKStar0')
        hSE = hSEMultVsKStar.ProjectionX('hSE')
        hSE.Rebin(bw)
        hSE.Sumw2()

        hMEMultVsKStar = inFile.Get(f'{comb}/ME/{region}/hCharmPtVsKStar0')
        hME = hMEMultVsKStar.ProjectionX('hME')
        hME.Rebin(bw)

        print(hSEMultVsKStar.GetEntries())
        print(hMEMultVsKStar.GetEntries())

        hCF = hSE/hME
        hCF.Scale(GetNormFactor(hSE, hME, 1500, 2000))
        hCF.SetName('hCF')

        hCF.Write()

print(f"output saved in {args.oFileName}")
