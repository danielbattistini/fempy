import argparse

from ROOT import TFile

import fempy

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFileCharm')
    parser.add_argument('inFileBeauty')
    parser.add_argument('oFile')
    parser.add_argument('--charmFraction', default=0.95)
    args = parser.parse_args()

    inFileCharm = TFile(args.inFileCharm)
    inFileBeauty = TFile(args.inFileBeauty)

    if set(fempy.utils.io.GetCombs(inFileCharm)) != set(fempy.utils.io.GetCombs(inFileBeauty)):
        fempy.error("different number of pairs")

    oFile = TFile(args.oFile, 'recreate')
    
    for comb in fempy.utils.io.GetCombs(inFileCharm):
        oFile.mkdir(comb)
        oFile.cd(comb)
        for event in ['SE', 'ME']:
            oFile.mkdir(f'{comb}/{event}')
            oFile.cd(f'{comb}/{event}')
            histsCharm = fempy.utils.io.GetHistsInDir(inFileCharm.Get(f'{comb}/{event}'))
            histsBeauty = fempy.utils.io.GetHistsInDir(inFileBeauty.Get(f'{comb}/{event}'))


            for hCharm, hBeauty in zip(histsCharm, histsBeauty):
                hQuark = hCharm.Clone()
                hQuark.Reset()
                hCharm.Sumw2()
                hBeauty.Sumw2()

                # hQuark.Add(hCharm, 0.5)
                hQuark.Add(hCharm, args.charmFraction)
                hQuark.Add(hBeauty, (1. - args.charmFraction)/hBeauty.GetEntries()*hCharm.GetEntries())
                hQuark.Write()

            regionsCharm = fempy.utils.GetRegions(inFileCharm.Get(f'{comb}/{event}'))
            regionsBeauty = fempy.utils.GetRegions(inFileBeauty.Get(f'{comb}/{event}'))
            print(regionsCharm, regionsBeauty)
            if set(regionsCharm) != set(regionsBeauty):
                fempy.error("not the same regions")

            for region in regionsCharm:
                oFile.mkdir(f'{comb}/{event}/{region}')
                oFile.cd(f'{comb}/{event}/{region}')
                histsCharm = fempy.utils.io.GetHistsInDir(inFileCharm.Get(f'{comb}/{event}/{region}'))
                histsBeauty = fempy.utils.io.GetHistsInDir(inFileBeauty.Get(f'{comb}/{event}/{region}'))
                print(histsCharm)
                for hCharm, hBeauty in zip(histsCharm, histsBeauty):
                    hQuark = hCharm.Clone()
                    hQuark.Reset()
                    hCharm.Sumw2()
                    hBeauty.Sumw2()

                    hQuark.Add(hCharm, args.charmFraction)
                    hQuark.Add(hBeauty, (1. - args.charmFraction)/hBeauty.GetEntries()*hCharm.GetEntries())
                    hQuark.Write()

    oFile.Close()
