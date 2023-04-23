import sys
import yaml

import argparse

from ROOT import TFile, TDatabasePDG

import fempy


def ComputeNormFactor(se, me, start, end):
    firstBin = se.FindBin(start*1.0001)
    lastBin = me.FindBin(end*0.9999)
    return me.Integral(firstBin, lastBin) / se.Integral(firstBin, lastBin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarPi', 'DstarK'))
    parser.add_argument('--syst', action='store_true', default=False)
    parser.add_argument('--stat', action='store_true', default=False)
    args = parser.parse_args()

    kStarBW = 50  # MeV/c
    if args.pair == 'DstarK':
        fitRanges = [[10, 450], [10, 400], [10, 500]]
        inFileData = TFile('/home/daniel/an/DstarK/2_luuksel/distr/Distr_data_nopc_kStarBW50MeV.root')
        inFileMC = TFile('~/an/DstarK/2_luuksel/distr/Distr_mchf_nopc_kStarBW50MeV_fromq.root')
        oFileName = f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_nopc_kStarBW50MeV_fromq.root'
        config = '/home/daniel/an/DstarK/cfg_gencf_DstarK_50MeV.yml'

        lpdg = 321

        weights1 = [0.78, 0.80, 0.77]
        radii1 = [0.86, 0.95, 0.79]
        radii2 = [2.03, 2.22, 1.91]
    else:
        print("not implemented")
        sys.exit()

    mLF = TDatabasePDG.Instance().GetParticle(lpdg).Mass()

    # load yaml file with lambda parameter
    with open(config, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    oFile = TFile(oFileName, 'recreate')
    for comb in ['sc', 'oc']:
        oFile.mkdir(comb)
        oFile.cd(comb)

        # Compute MC CF
        hSEMC = inFileMC.Get(f'{comb}/SE/sgn/hCharmMassVsKStar0').ProjectionX('hSEMC')
        hMEMC = inFileMC.Get(f'{comb}/ME/sgn/hCharmMassVsKStar0').ProjectionX('hMEMC')

        rebinFactor = round(kStarBW/hSEMC.GetBinWidth(1))
        hSEMC.Rebin(rebinFactor)
        hMEMC.Rebin(rebinFactor)
        hSEMC.Scale(ComputeNormFactor(hSEMC, hMEMC, 1000, 1500))
        hCFMC = hSEMC/hMEMC
        hSEMC.Write()
        hMEMC.Write()
        hCFMC.SetName('hCFMC')
        hCFMC.Write()

        regions = fempy.utils.GetRegions(inFileData.Get('sc/SE'))
        dCFData = {}
        for region in regions:
            hSEData = inFileData.Get(f'{comb}/SE/{region}/hCharmMassVsKStar0').ProjectionX(f'hSE_{region}')
            hMEData = inFileData.Get(f'{comb}/ME/{region}/hCharmMassVsKStar0').ProjectionX(f'hME_{region}')

            rebinFactor = round(kStarBW/hSEData.GetBinWidth(1))
            hSEData.Rebin(rebinFactor)
            hMEData.Rebin(rebinFactor)
            hSEData.Scale(ComputeNormFactor(hSEData, hMEData, 1000, 1500))
            hCFData = hSEData/hMEData
            hSEData.Write()
            hMEData.Write()
            hCFData.SetName(f'hCF_{region}')
            hCFData.Write()
            dCFData[region] = hCFData

    print(f'output saved in {oFileName}')
    oFile.Close()
