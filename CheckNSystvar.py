import argparse

from ROOT import TFile, TCanvas, TH1D, gROOT

import fempy


def GetNVars(inFile):
    vars = []
    for name in fempy.utils.io.GetObjNamesInDir(inFile.Get('pp/SE')):
        if name.find('_') >= 0:
            vars.append(name[:name.find('_')])
    return len(list(dict.fromkeys(vars)))


def CheckNSystVar(inFile, oDir, suffix, pair):
    inFile = TFile.Open(inFile)
    for comb in ['pp', 'mm', 'pm', 'mp', 'sc', 'oc']:
        femtoEntriesList = []
        nVars = GetNVars(inFile)
        for syst in range(nVars):
            hTmp = inFile.Get(f'{comb}/SE/hCharmMassVsKStar{syst}').ProjectionX()
            femtoEntriesList.append(hTmp.Integral(1, hTmp.GetXaxis().FindBin(200*0.9999)))
            print(hTmp.GetXaxis().FindBin(200*0.999))
        print(femtoEntriesList)

        gROOT.SetBatch(False)
        cVars = TCanvas(f'cVars_{comb}_SE', f'cVars_{comb}_SE', 600, 600)
        title = fempy.utils.format.TranslateToLatex(f'k{pair}_{comb};Variation;Entries for k*<200 MeV/c')
        hVars = TH1D(f'hVars_{comb}_SE', title, nVars, -0.5, nVars-0.5)
        for syst, femtoEntries in enumerate(femtoEntriesList):
            hVars.SetBinContent(syst+1, femtoEntries)
        hVars.Draw()

        cVars.SetLeftMargin(0.15)
        cVars.SaveAs(f'{oDir}/cVars_{comb}_SE.pdf')
        cVars.SaveAs(f'{oDir}/cVars_{comb}_SE.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('oDir')
    parser.add_argument('--suffix', default='')
    parser.add_argument('--pair', default='')
    args = parser.parse_args()

    CheckNSystVar(args.inFile, args.oDir, args.suffix, pair=args.pair)
