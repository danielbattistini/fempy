# https://link.springer.com/content/pdf/10.1007/JHEP01(2022)106.pdf

import sys
import argparse
import os

from ROOT import TFile, TH1, TH1D, TCanvas, TLatex, TLegend, kRed, gStyle, kBlue, kRed

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

oFileName = os.path.join(oDir, "MultTuning_pp13TeV.root")
oFile = TFile(oFileName, "recreate")

targetAvgMult = 31.5 # taken from https://link.springer.com/content/pdf/10.1007/JHEP01(2022)106.pdf Tab. 2

for iMult in dMult:
    dMult[iMult].SetLineWidth(2)
    dMult[iMult].Write()
    
    if (dMult[iMult].GetMean() > targetAvgMult):
        print(f"target mult matched for V0Mult emul > {iMult}:", dMult[iMult].GetMean())

oFile.Close()

