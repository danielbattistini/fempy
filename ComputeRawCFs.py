import argparse
import os

from ROOT import TFile

import fempy
from fempy import CorrelationFunction
from fempy.utils import Pair

parser = argparse.ArgumentParser()
parser.add_argument('--pair', default="")
parser.add_argument('inFileName')
parser.add_argument('oDir')
parser.add_argument('--suffix', default="")
args = parser.parse_args()

inFile = TFile(args.inFileName)
oFileName = os.path.join(args.oDir, "RawCF.root" if args.suffix == '' else f'RawCF_{args.suffix}.root')
oFile = TFile(oFileName, 'recreate')

combs = fempy.utils.io.GetSubdirsInDir(inFile)
regions = fempy.utils.io.GetSubdirsInDir(inFile.Get(f'{combs[0]}/SE'))
pair = Pair(args.pair)

binWidths = [1, 2, 4, 5, 10, 15, 20, 40, 50]  # MeV/c

for comb in combs:
    for region in regions:
        hSEMultVsKStar = inFile.Get(f'{comb}/SE/{region}/hMultVsKStar0')
        hMEMultVsKStar = inFile.Get(f'{comb}/ME/{region}/hMultVsKStar0')

        print(hSEMultVsKStar.GetEntries())
        print(hMEMultVsKStar.GetEntries())

        hSEProj = hSEMultVsKStar.ProjectionX()
        hMEProj = hMEMultVsKStar.ProjectionX()

        for bw in binWidths:
            hSE = hSEProj.Clone()
            hSE.Rebin(bw)

            hME = hMEProj.Clone()
            hME.Rebin(bw)
            print(hMEMultVsKStar)

            print(pair.norm_range)
            hCF = CorrelationFunction(se=hSE, me=hME, norm=pair.norm_range).get_cf()
            hCF.Write()

print(f"output saved in {oFileName}")
