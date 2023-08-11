# https://link.springer.com/content/pdf/10.1007/JHEP01(2022)106.pdf

import sys
import argparse
import os

from ROOT import TFile, gStyle

gStyle.SetPalette(55)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileName', metavar='text')
parser.add_argument('oFileName', metavar='text')
args = parser.parse_args()

inFile = TFile.Open(args.inFileName)
oFileName = args.oFileName

if inFile == None:
    print("Error: the file could not be opened. Exit!")
    sys.exit()
    
hMult = inFile.Get("hMult")
nBinsV0 = hMult.GetNbinsX()

dMult = {}
for iV0Bin in range(nBinsV0):
    iMult = int(hMult.GetXaxis().GetBinCenter(iV0Bin))
    dMult[iMult] = hMult.ProjectionY(f"hMult_{iMult}", iMult, -1) # integrate until last bin    
    dMult[iMult].Scale(1./dMult[iMult].GetEntries())

oFile = TFile(oFileName, "recreate")

targetAvgMult = 31.5 # taken from https://link.springer.com/content/pdf/10.1007/JHEP01(2022)106.pdf Tab. 2

for iMult, hMult in dMult.items():
    hMult.SetLineWidth(2)
    hMult.Write()

    if hMult.GetMean() > targetAvgMult:
        print(f"target mult matched for V0Mult emul > {iMult}:", hMult.GetMean())

oFile.Close()
