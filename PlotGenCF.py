import argparse

from ROOT import (
    TFile,
    TCanvas,
    TLatex,
    gStyle,
    kBlue,
    kGray,
    kOrange,
    TLegend,
    TGraphErrors,
    kWhite,
    kRed,
    kGreen,
    kBlack
)

from fempy.utils.format import TranslateToLatex


def LoadObjects(pair, suffix):
    objs = {
        'sc': {},
        'oc': {},
    }

    if pair == 'DstarK':
        file = TFile(f'/home/daniel/an/DstarK/2_luuksel/GenCFCorr_{suffix}.root')
        for comb in ['sc', 'oc']:
            objs[comb] = {
                'hstat': file.Get(f'{comb}/hCFGenStat'),
                'stat': file.Get(f'{comb}/gCFGenStat'),
                'syst': file.Get(f'{comb}/gCFGenSyst'),
                'llstat': file.Get(f'{comb}/gLLStat'),
                'lltot': file.Get(f'{comb}/gLLTot'),
                'coulomb': file.Get(f'{comb}/gCoulomb'),
            }
            objs[comb]['hstat'].SetDirectory(0)
        objs['sc']['a0'] = -0.305
        objs['sc']['a0stat'] = 0.235
        objs['sc']['a0syst'] = 0.083
        objs['oc']['a0'] = -0.232
        objs['oc']['a0stat'] = 0.181
        objs['oc']['a0syst'] = 0.062
        config = {
            'pltRangeX': [0, 299.9999],
            'pltRangeY': [-0.2, 3.4999],
        }
        file.Close()
    elif pair == 'DstarPi':
        file = TFile(f'/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_{suffix}.root')
        for comb in ['sc', 'oc']:
            objs[comb] = {
                'hstat': file.Get(f'{comb}/hCFGenStat'),
                'stat': file.Get(f'{comb}/gCFGenStat'),
                'syst': file.Get(f'{comb}/gCFGenSyst'),
                'llstat': file.Get(f'{comb}/gLLStat'),
                'lltot': file.Get(f'{comb}/gLLTot'),
                'coulomb': file.Get(f'{comb}/gCoulomb'),
            }
            objs[comb]['hstat'].SetDirectory(0)

        # old values obtained from non simul fit
        objs['sc']['a0'] = 0.112
        objs['sc']['a0stat'] = 0.035
        objs['sc']['a0syst'] = 0.026
        objs['oc']['a0'] = -0.047
        objs['oc']['a0stat'] = 0.033
        objs['oc']['a0syst'] = 0.019

        # from simul fit
        objs['3/2'] = {}
        objs['3/2']['a0'] = 0.12
        objs['3/2']['a0stat'] = 0.04
        objs['3/2']['a0syst'] = 0.02
        objs['1/2'] = {}
        objs['1/2']['a0'] = -0.12
        objs['1/2']['a0stat'] = 0.06
        objs['1/2']['a0syst'] = 0.03

        config = {
            'pltRangeX': [0, 300],
            'pltRangeY': [0.7, 1.9],
        }
        file.Close()

    return objs, config


def Setstyle():
    gStyle.SetTextFont(42)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)

    gStyle.Reset("Plain")
    gStyle.SetOptTitle(False)
    gStyle.SetTitleBorderSize(0)
    gStyle.SetOptStat(0)
    gStyle.SetCanvasColor(10)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetFrameLineWidth(1)
    gStyle.SetFrameFillColor(kWhite)
    gStyle.SetPadColor(10)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetHistLineWidth(1)
    gStyle.SetHistLineColor(kRed)
    gStyle.SetFuncWidth(2)
    gStyle.SetFuncColor(kGreen)
    gStyle.SetLineWidth(2)
    gStyle.SetLabelSize(0.045, "xyz")
    gStyle.SetLabelOffset(0.01, "y")
    gStyle.SetLabelOffset(0.01, "x")
    gStyle.SetLabelColor(kBlack, "xyz")
    gStyle.SetTitleSize(0.05, "xyz")
    gStyle.SetTitleOffset(1.25, "y")
    gStyle.SetTitleOffset(1.2, "x")
    gStyle.SetTitleFillColor(kWhite)
    gStyle.SetTextSizePixels(26)
    gStyle.SetTextFont(42)
    gStyle.SetLegendFillColor(kWhite)
    gStyle.SetLegendFont(42)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetErrorX(0.005)
    # gStyle.SetPalette(kBird)


def SetHistStyle(hist):
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelOffset(0.01)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelOffset(0.01)
    hist.GetYaxis().SetTitleOffset(1.25)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)


def ComputeBinBrackets(hist):
    gBrackets = TGraphErrors(1)
    print(hist)
    for iBin in range(hist.GetNbinsX()):
        gBrackets.SetPoint(iBin, hist.GetBinCenter(iBin+1), hist.GetBinContent(iBin+1))
        gBrackets.SetPointError(iBin, hist.GetBinWidth(iBin+1)/2, 0)
    return gBrackets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', choices=('DstarK', 'DstarPi'))
    parser.add_argument('--simfit', default=False, action='store_true')
    args = parser.parse_args()

    if args.pair == 'DstarPi':
        suffix = 'nopc_kStarBW50MeV_bs5000syst'
    else:
        suffix = 'nopc_kStarBW50MeV_fromq_bs5000syst'

    Setstyle()

    objects, cfg = LoadObjects(args.pair, suffix)

    cCF = TCanvas('cCF', '', 1000, 500)

    cCF.Divide(2, 1)
    for iPad, comb in enumerate(['sc', 'oc']):
        objs = objects[comb]
        pad = cCF.cd(iPad+1)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.03)
        pad.SetFrameLineWidth(1)
        pad.DrawFrame(
            cfg['pltRangeX'][0],
            cfg['pltRangeY'][0],
            cfg['pltRangeX'][1],
            cfg['pltRangeY'][1], TranslateToLatex(';__kStarMeV__;__C__'))

        objs['lltot'].SetFillColorAlpha(kBlue+1, 0.7)
        objs['lltot'].SetLineColorAlpha(0, 0)
        objs['lltot'].Draw('same e3')
        objs['fit'] = objs['lltot'].Clone()
        objs['fit'].SetName('gFit')
        objs['fit'].Write('gFit')

        # objs['llstat'].SetFillColorAlpha(kBlue+1, 0.7)
        # objs['llstat'].SetLineColorAlpha(0, 0)
        # objs['llstat'].Draw('same e3')

        objs['coulomb'].SetLineColor(kOrange+7)
        objs['coulomb'].SetLineStyle(1)
        objs['coulomb'].SetLineWidth(2)
        objs['coulomb'].Draw('same l')

        objs['syst'].SetFillColorAlpha(kGray+2, 0.65)
        objs['syst'].Draw('same e2')

        objs['stat'].SetMarkerSize(1)
        objs['stat'].SetMarkerStyle(24)
        objs['stat'].Draw('same pez')

        gBrackets = ComputeBinBrackets(objs['hstat'])
        SetHistStyle(gBrackets)
        gBrackets.SetLineWidth(1)
        gBrackets.SetLineColorAlpha(kBlack, 0.9)
        gBrackets.DrawClone("same []")

        tl = TLatex()
        tl.SetTextSize(0.04)
        tl.SetNDC(True)
        if iPad == 0:
            tl.DrawLatex(0.22, 0.88, 'ALICE pp #sqrt{#it{s}} = 13 TeV')
            tl.DrawLatex(0.22, 0.82, 'High-mult. (0 #minus 0.17% INEL > 0)')
            tl.DrawLatex(0.22, 0.76, TranslateToLatex(f'k{args.pair}_{comb}'))

            leg = TLegend(0.2, 0.55, 0.9, 0.7)
            leg.SetTextSizePixels(12)
            leg.SetTextSize(0.035)
            leg.SetLineColorAlpha(0, 0)
            leg.SetFillStyle(0)
            leg.AddEntry(objs['stat'], 'Genuine CF', 'pel')
            leg.AddEntry(objs['coulomb'], 'Coulomb only', 'l')
            leg.AddEntry(objs['lltot'], 'Lednický-Lyuboshits model', 'f')
            leg.Draw()

            if not args.simfit:
                text = f'a_{{0}} = {objs["a0"]:.2f} #pm {objs["a0stat"]:.2f} (stat) #pm {objs["a0syst"]:.2f} (syst) fm'
                tl.DrawLatex(0.20, 0.2, text)
        else:
            tl.DrawLatex(0.20, 0.88, TranslateToLatex(f'k{args.pair}_{comb}'))

            if args.simfit:
                text3_2 = (
                    f'a_{{0}}(I=3/2) = {objects["3/2"]["a0"]:.2f} '
                    f' #pm {objects["3/2"]["a0stat"]:.2f} (stat)'
                    f' #pm {objects["3/2"]["a0syst"]:.2f} (syst) fm'
                )
                tl.DrawLatex(0.20, 0.25, text3_2)

                text1_2 = (
                    f'a_{{0}}(I=1/2) = {objects["1/2"]["a0"]:.2f}'
                    f' #pm {objects["1/2"]["a0stat"]:.2f} (stat)'
                    f' #pm {objects["1/2"]["a0syst"]:.2f} (syst) fm'
                )
                tl.DrawLatex(0.20, 0.25-0.05, text1_2)
            else:
                text = f'a_{{0}} = {objs["a0"]:.2f} #pm {objs["a0stat"]:.2f} (stat) #pm {objs["a0syst"]:.2f} (syst) fm'
                tl.DrawLatex(0.20, 0.2, text)

    for ext in ['pdf', 'png', 'eps', 'root']:
        if args.pair == 'DstarK':
            cCF.SaveAs(f'/home/daniel/an/DstarK/2_luuksel/GenCFCorrPlot_{args.pair}_{suffix}.{ext}')
        else:
            oFileName = (
                '/home/daniel/an/DstarPi/20_luuksel/'
                f'GenCFCorrPlot_{args.pair}_{suffix}{"_simfit" if args.simfit else ""}.{ext}'
            )
            cCF.SaveAs(oFileName)


if __name__ == '__main__':
    main()
