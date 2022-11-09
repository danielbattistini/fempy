'''Compute the raw correlation functions with the output of FemtoDream'''

import sys
import os

import argparse

from ROOT import TFile, TCanvas, gROOT, TH1I

from core.CorrelationFunction import CorrelationFunction
from utils.format import TranslateToLatex, GetNormRangeFromLabel
from utils.io import GetObjectFromFile


def GetListOfPairKeys(file):
    '''Returns the list of types of CF as in FD: signal, sideband, MCTruth etc'''
    keys = [key.GetName() for key in list(file.GetListOfKeys())]
    keys = [key.replace('HM_CharmFemto_', '').replace("_Results0", '') for key in keys if 'Results' in key]

    return keys


def GetPairLabel(key):
    '''Converts FD labels for particles into fempy standard'''
    if 'DstarPi' in key:
        return 'DstarPi'
    elif 'Dkaon' in key:
        print("ciccio")
        return 'DK'
    elif 'Dpion' in key or 'DplusPion' in key:
        print("pasticcio")
        return 'DPi'
    else:
        print("Error: pair not implemented. Exit!")
        sys.exit()


def GetPairRegion(key):
    '''Convert FD labels for mass regions to fempy standard'''
    if 'DstarPion' == key or 'Dkaon' == key or 'Dpion' == key or 'DplusPion' == key:
        return 'sgn'
    elif 'SBRight_' in key:
        return 'sbr'
    elif 'SBLeft_' in key:
        return 'sbl'
    else:
        print("Error: region not implemented. Exit!")
        sys.exit()


gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName')
parser.add_argument('oFileName')
args = parser.parse_args()

inFile = TFile(args.inFileName)
pairsKeys = GetListOfPairKeys(inFile)
print(pairsKeys)

oFile = TFile(os.path.splitext(args.oFileName)[0]+'.root', 'recreate')

binWidths = [5, 10, 20, 50]  # MeV
dSE, dME, dCF = ({} for _ in range(3))

# load ME
cCF = TCanvas("cCF", "cCF", 600, 600)
cCF.SaveAs(os.path.splitext(args.oFileName)[0] + '.pdf')

for key in pairsKeys:
    basename = f"HM_CharmFemto_{key}_Results0/" * 2

    pairLabel = GetPairLabel(key)
    pairRegion = GetPairRegion(key)
    print(key, pairLabel, pairRegion)

    # load SE
    objName = basename + "Particle0_Particle2/SEDist_Particle0_Particle2"
    dSE[f"hSE_{pairLabel}_pp_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle1_Particle3/SEDist_Particle1_Particle3"
    dSE[f"hSE_{pairLabel}_mm_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle0_Particle3/SEDist_Particle0_Particle3"
    dSE[f"hSE_{pairLabel}_mp_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle1_Particle2/SEDist_Particle1_Particle2"
    dSE[f"hSE_{pairLabel}_pm_{pairRegion}"] = GetObjectFromFile(inFile, objName)

    # load ME
    objName = basename + "Particle0_Particle2/MEDist_Particle0_Particle2"
    dME[f"hME_{pairLabel}_pp_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle1_Particle3/MEDist_Particle1_Particle3"
    dME[f"hME_{pairLabel}_mm_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle0_Particle3/MEDist_Particle0_Particle3"
    dME[f"hME_{pairLabel}_mp_{pairRegion}"] = GetObjectFromFile(inFile, objName)
    objName = basename + "Particle1_Particle2/MEDist_Particle1_Particle2"
    dME[f"hME_{pairLabel}_pm_{pairRegion}"] = GetObjectFromFile(inFile, objName)

    # sum same charge
    keySESum = f"hSE_{pairLabel}_sc_{pairRegion}"
    dSE[keySESum] = dSE[f"hSE_{pairLabel}_pp_{pairRegion}"].Clone()
    dSE[keySESum].Add(dSE[f"hSE_{pairLabel}_mm_{pairRegion}"])
    keyMESum = f"hME_{pairLabel}_sc_{pairRegion}"
    dME[keyMESum] = dME[f"hME_{pairLabel}_pp_{pairRegion}"].Clone()
    dME[keyMESum].Add(dME[f"hME_{pairLabel}_mm_{pairRegion}"])

    # sum opposite charge
    keySESum = f"hSE_{pairLabel}_oc_{pairRegion}"
    dSE[keySESum] = dSE[f"hSE_{pairLabel}_pm_{pairRegion}"].Clone()
    dSE[keySESum].Add(dSE[f"hSE_{pairLabel}_mp_{pairRegion}"])
    keyMESum = f"hME_{pairLabel}_oc_{pairRegion}"
    dME[keyMESum] = dME[f"hME_{pairLabel}_pm_{pairRegion}"].Clone()
    dME[keyMESum].Add(dME[f"hME_{pairLabel}_mp_{pairRegion}"])

    hNPairs = TH1I(f"hNPairs_kstarlt200MeV_{pairLabel}_{pairRegion}", 'hNPairs_kstarlt200MeV', 6, -0.5, 5.5)
    for iComb, comb in enumerate(['pp', 'mm', 'sc', 'pm', 'mp', 'oc']):
        se = dSE[f"hSE_{pairLabel}_{comb}_{pairRegion}"]
        nPairs = se.Integral(1, se.FindBin(0.2-1.e-6))
        hNPairs.SetBinContent(iComb + 1, nPairs)
        hNPairs.GetXaxis().SetBinLabel(iComb + 1, TranslateToLatex(f'k{pairLabel}_{comb}'))

    hNPairs.Write()

    # compute CF
    binWidth = dME[keyMESum].GetXaxis().GetBinWidth(1)

    # write to file
    for bwNew in binWidths:
        for keySE, keyME in zip(dSE, dME):
            if 'DPi' in keySE:
                norm = GetNormRangeFromLabel('DPi')
            elif 'DK' in keySE:
                norm = GetNormRangeFromLabel('DK')
            elif 'DstarPi' in keySE:
                norm = GetNormRangeFromLabel('DstarPi')
            else:
                print("Error: pair not implemented. Exit!")
                sys.exit()

            cf = CorrelationFunction(dSE[keySE], dME[keyME], um='GeV2MeV', norm=norm)
            cf.Rebin(int(bwNew/1000/binWidth))
            dCF[f"hCF_{keySE[4:]}_bw{bwNew}MeV"] = cf.GetCF()

    for keyCF in dCF:
        dCF[keyCF].SetTitle(keyCF)
        dCF[keyCF].Draw("pe")

        cCF.SaveAs(os.path.splitext(args.oFileName)[0]+'.pdf')

    for cf in dCF:
        dCF[cf].Write(cf)

cCF.SaveAs(os.path.splitext(args.oFileName)[0] + '.pdf')
oFile.Close()
print(f"Output saved to {os.path.splitext(args.oFileName)[0]}.root")
