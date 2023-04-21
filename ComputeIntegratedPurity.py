import argparse

import numpy as np

import fempy

from ROOT import TFile, TCanvas, gInterpreter, gPad, TLatex
gInterpreter.ProcessLine('#include "fempy/MassFitter.hxx"')
from ROOT import MassFitter

def ComputeIntegratedPurity(inFileName, suffix):
    inFile = TFile.Open(inFileName)

    suffix = f'_{suffix}' if suffix != '' else suffix
    oFileName = f'Purity{suffix}.root'
    oFile = TFile(oFileName, 'recreate')


    fitRange = [0.141, 0.154]

    for comb in ['sc', 'oc']:
        oFile.mkdir(comb)
        oFile.cd(comb)

        hMass = inFile.Get(f'{comb}/SE/hCharmMassVsKStar0').ProjectionY()
        firstBin = hMass.GetXaxis().FindBin(fitRange[0]*1.0001)
        lastBin = hMass.GetXaxis().FindBin(fitRange[1]*0.9999)

        print(firstBin, lastBin)

        particleYield = hMass.Integral(firstBin, lastBin)
        fitter = MassFitter(hMass, 'gaus', 'powex', fitRange[0], fitRange[1])

        fitter.Fit()

        # chi2 = fitter.GetChi2Ndf()
        sgn = fitter.GetSignal(2, 'data_minus_bkg')
        sgnUnc = fitter.GetSignalUnc(2, 'data_minus_bkg')
        
        bkg = fitter.GetBackground(2)
        bkgUnc = fitter.GetBackgroundUnc(2)

        soverb = sgn/bkg
        soverbUnc = soverb * np.sqrt((sgnUnc/sgn)**2 + (bkgUnc/bkg)**2)

        purity = sgn / (sgn + bkg)
        purityUnc = purity * np.sqrt((sgnUnc/sgn)**2 + (1./np.sqrt(particleYield))**2)

        
        # draw
        cPurity = TCanvas(f"cPurity_{comb}", "Purity {comb}", 600, 600)
        fitter.Draw(gPad, 'data_minus_bkg')

        tl = TLatex()
        tl.SetNDC()
        tl.SetTextSize(0.035)
        tl.SetTextFont(42)
        tl.DrawLatex(0.6, 0.85, f'p={purity:.2f} #pm {purityUnc:.2f}')
        tl.DrawLatex(0.6, 0.85 - 0.05, f'S/B={soverb:.2f} #pm {soverbUnc:.2f}')

        
        cPurity.SaveAs('test.png')
        cPurity.Write()
    print(f"output saved in {oFileName}")
    oFile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('--suffix', default='')
    args = parser.parse_args()
    
    ComputeIntegratedPurity(args.inFile, args.suffix)


