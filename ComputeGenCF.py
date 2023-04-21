import yaml
import sys
import argparse

import fempy
from fempy import CorrelationFunction

from ROOT import TFile, gROOT


def IsLamParMatValid(lam_par):
    total = 0
    for _, h_lam in lam_par.items():
        for _, l_lam in h_lam.items():
            total += l_lam
    return abs(total - 1) < 1e-6


def SumLamPar(lam_par, treamtments):
    lam_par_summed = {
        'gen': 0,
        'flat': 0,
    }
    for hk, h_lam in lam_par.items():
        for lk, l_lam in h_lam.items():
            lam_par_summed[treamtments[hk][lk]] += l_lam
    return lam_par_summed


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('cfg')
    parser.add_argument('oFile')
    parser.add_argument('--dry', default=False, action='store_true')
    args = parser.parse_args()

    with open(args.cfg, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    lam_par_matr = {}
    for h_contrib in cfg['heavy']:
        h_key = list(h_contrib.keys())[0]
        h_purity = h_contrib[h_key]['purity']
        h_frac = h_contrib[h_key]['frac']
        lam_par_matr[h_key] = {}
        for l_contrib in cfg['light']:
            l_key = list(l_contrib.keys())[0]
            l_purity = l_contrib[l_key]['purity']
            l_frac = l_contrib[l_key]['frac']

            lam_par_matr[h_key][l_key] = l_frac * l_purity * h_frac * h_purity
    print('\n\n\n---> ', lam_par_matr, 'sum=', IsLamParMatValid(lam_par_matr))
    lam_par_sum = SumLamPar(lam_par_matr, cfg['treatment'])
    if args.dry:
        print(lam_par_sum)
        sys.exit()

    inFile = TFile(args.inFile)
    if not args.dry:
        oFile = TFile(args.oFile, "recreate")
    # inFile = TFile("~/an/DstarPi/18_fixmix/cf/RawCF_data_nopc_kStarBW15MeV_gaus_charmMassBW0.2MeV_lowKStarCharmMassBW0.4MeV_until50MeV.root")
    # oFileName = "~/an/DstarPi/18_fixmix/cf/GenCF_nopc_kStarBW15MeV_gaus_charmMassBW0.2MeV_lowKStarCharmMassBW0.4MeV_until50MeV.root"

    histNames = fempy.utils.io.GetHistNamesInDir(inFile.Get("pp"))
    nSysts = len([name for name in histNames if "hCFPurityRewFromDataMinusBkg" in name])
    if nSysts == 0:
        nSysts = len([name for name in histNames if "hCFstd_sgn" in name])

    print(nSysts)
    # gROOT.SetBatch(False)
    # print(systs)

    
    for comb in ['sc', 'oc']:
        if not args.dry:
            oFile.mkdir(comb)

        CFMinijet = 1
        if cfg['mj']['enable']:
            inFileMJ = TFile.Open(cfg['mj']['file'])
            CFMinijet = inFileMJ.Get(f'{comb}/hCFstd_sgn0')

        for syst in range(nSysts):
            hCFRaw = inFile.Get(f"{comb}/hCFPurityRewFromDataMinusBkg{syst}")
            if not hCFRaw:
                hCFRaw = inFile.Get(f"{comb}/hCFstd_sgn{syst}")
            if not hCFRaw:
                print("hCFRaw not leaded")
                inFile.ls()
                fempy.error("hCFRaw not leaded")
            # print(hCFRaw.GetNbinsX())
            CFRaw = CorrelationFunction(cf=hCFRaw)
            # print("kakas")
            print(CFRaw, CFMinijet,  lam_par_sum['flat'], lam_par_sum['gen'])
            # CFMinijet.Draw()
            CFGen = (CFRaw/CFMinijet - lam_par_sum['flat'])/lam_par_sum['gen']
            # CFGen = (CFRaw - lam_par_sum['flat'])/lam_par_sum['gen']/CFMinijet

            CFGen.get_cf().Draw()
            hCFRaw.SetLineColor(3)
            hCFRaw.Draw('same')

            if not args.dry:
                oFile.cd(comb)
                CFGen.get_cf().Write(f"hGenCF{ſyst}")

                # save graph
                if syst == 0:
                    hCF = CFGen.get_cf()
                    gCF = inFile.Get(f"{comb}/gCFPurityRewFromDataMinusBkg")
                    inFile.Get(comb).ls()
                    if gCF == None:
                        gCF = inFile.Get(f"{comb}/gCFstd_sgn")

                    yCF = [hCF.GetBinContent(iBin+1) for iBin in range(hCF.GetNbinsX())]
                    yCFUnc = [hCF.GetBinError(iBin+1) for iBin in range(hCF.GetNbinsX())]
                    xCF = [gCF.GetPointX(iBin) for iBin in range(hCF.GetNbinsX())]
                    xCFUnc = [gCF.GetErrorX(iBin) for iBin in range(hCF.GetNbinsX())]
                    
                    print(xCF)
                    print(yCF)
                    print(xCFUnc)
                    print(yCFUnc)
                    for iPoint, (x, y, xUnc, yUnc) in enumerate(zip(xCF, yCF, xCFUnc, yCFUnc)):
                        gCF.SetPoint(iPoint, x, y)
                        gCF.SetPointError(iPoint, xUnc, yUnc)
                    gCF.Write(f"gGenCF{ſyst}")

    if not args.dry:
        
        print("output saved in ", args.oFile)
        oFile.Close()
