import sys
import numpy as np
import sys
import os
import argparse

from rich import print

from ROOT import TFile, TCanvas, kRed, TLegend, gStyle, gPad, TLatex, SetOwnership, kRed, kOrange, kMagenta, kViolet, kAzure, gROOT

sys.path.append('../')
from fempy.utils.format import GetNormRangeFromPDG, GetNormRangeFromPDG, dPDG2Label, TranslateToLatex, FigInit
from correlation_function import CorrelationFunction

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName')
parser.add_argument('oFileName')
args = parser.parse_args()

FigInit()

# dCombo2Latex = {
#     'DplusPi': {
#         'pp': 'D^{+}#pi^{+}',
#         'mm': 'D^{#minus}#pi^{#minus}',
#         'sc': 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}',
#         'oc': 'D^{+}#pi^{#minus} #oplus D^{#minus}#pi^{+}'
#     },
#     'DplusK': {
#         'pp': 'D^{+}K^{+}',
#         'mm': 'D^{#minus}K^{#minus}',
#         'sc': 'D^{+}K^{+} #oplus D^{#minus}K^{#minus}',
#         'oc': 'D^{+}K^{#minus} #oplus D^{#minus}K^{+}'
#     },
#     'DstarPi': {
#         'pp': 'D*^{+}#pi^{+}',
#         'mm': 'D*^{#minus}#pi^{#minus}',
#         'sc': 'D*^{+}#pi^{+} #oplus D*^{#minus}#pi^{#minus}',
#         'oc': 'D*^{+}#pi^{#minus} #oplus D*^{#minus}#pi^{+}'
#     },
#     'DstarK': {
#         'pp': 'D*^{+}K^{+}',
#         'mm': 'D*^{#minus}K^{#minus}',
#         'sc': 'D*^{+}K^{+} #oplus D*^{#minus}K^{#minus}',
#         'oc': 'D*^{+}K^{#minus} #oplus D*^{#minus}K^{+}'
#     },
# }


# set style
gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLineWidth(2)

gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetPadBottomMargin(0.1)


rebin = 50
lightPDG = [211, 321]  # pions, kaons
charmPDG = [411]  # D+, D*
suffix = ''

charmOrigins = ['prompt', 'nonprompt', 'any']

# load SE/ME
dSE = {}
dME = {}

c2D = TCanvas("c2D", "c2D", 1800, 800)
c2D.Divide(3, 1)
# load SE/ME, compute CF
# inFileName = '/home/daniel/an/DKDpi/mj/cfvspt/prod/pythia_HMpp13TeV_tune-kMonash_process-kSoftQCD_distr.root'
# inFileName = '/home/daniel/an/DKDpi/mj/cfvspt/5_promptnonpromptsep/pythiamj_tune-Monash_process-SoftQCD.root'
inFile = TFile(args.inFileName)
oFileNameBase = os.path.splitext(args.oFileName)[0]

for cPDG in charmPDG:
    for lPDG in lightPDG:
        for part in ['', '-']:
            comb = 'sc' if part == '' else 'oc'
            for charmOrigin in ['prompt', 'nonprompt']:

                se = inFile.Get(f'hPairSE_{cPDG}_{part}{lPDG}_{charmOrigin}')
                print(se)
                if se == None:
                    print(f'The histogram \'hPairSE_{cPDG}_{part}{lPDG}_{charmOrigin}\' couldn\'t be loaded. Exit!')
                    sys.exit()

                se.SetDirectory(0)
                dSE[f'hSE_{cPDG}_{part}{lPDG}_{charmOrigin}'] = se

                me = inFile.Get(f'hPairME_{cPDG}_{part}{lPDG}_{charmOrigin}')
                me.SetDirectory(0)
                dME[f'hME_{cPDG}_{part}{lPDG}_{charmOrigin}'] = me

                pad = c2D.cd(1)
                pad.SetRightMargin(0.15)
                pad.SetTopMargin(0.1)
                se.SetTitle("same-event")
                se.Draw('colz')

                # c2D.cd()
                # title = TLatex()
                # title.SetTextAlign(13)
                # title.SetNDC()
                # title.SetTextSize(0.05)
                # title.SetTextFont(42)
                # title.DrawLatex(0.4, 0.95, f'{dCombo2Latex[dPDG2Label[cPDG]}{dPDG2Label[lPDG]][comb]} {charmOrigin}')

                pad = c2D.cd(2)
                pad.SetTopMargin(0.1)
                pad.SetRightMargin(0.15)
                me.SetTitle("mixed-event")
                me.Draw('colz')

                pad = c2D.cd(3)
                pad.SetTopMargin(0.1)
                pad.SetRightMargin(0.15)
                cf = se.Clone(f'hCF_{cPDG}_{part}{lPDG}_{charmOrigin}')
                cf.Scale(1./cf.GetEntries())
                cf.RebinX(50)
                me2 = me.Clone()
                me2.Scale(1./me2.GetEntries())
                me2.RebinX(50)

                cf.Divide(me2)
                cf.GetZaxis().SetTitle('#it{C}(#it{k}*)')
                cf.SetTitle(TranslateToLatex(f'CF k{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb} {charmOrigin}'))
                # cf.SetTitle(f"CF {dCombo2Latex[dPDG2Label[cPDG]}{dPDG2Label[lPDG]][comb]} {charmOrigin}")
                cf.Draw('colz')
                c2D.SaveAs(f'{oFileNameBase}_2dcf_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}.png')
                # sys.exit()

inFile.Close()

# project CFs in momentum
ptMins = [1]
ptMaxs = [10]

# include the integrated ptbin
ptMins.append(ptMins[0])
ptMaxs.append(ptMaxs[-1])

# sum prompt and nonprompt
for cPDG in charmPDG:
    for lPDG in lightPDG:
        for part in ['', '-']:
            comb = 'sc' if part == '' else 'oc'

            # SE
            prompt = dSE[f'hSE_{cPDG}_{part}{lPDG}_prompt']
            nonprompt = dSE[f'hSE_{cPDG}_{part}{lPDG}_nonprompt']

            any = prompt.Clone(f'hSE_{cPDG}_{part}{lPDG}_any')
            any.Sumw2()
            any.Add(nonprompt)
            dSE[f'hSE_{cPDG}_{part}{lPDG}_any'] = any

            # ME
            prompt = dME[f'hME_{cPDG}_{part}{lPDG}_prompt']
            nonprompt = dME[f'hME_{cPDG}_{part}{lPDG}_nonprompt']

            any = prompt.Clone(f'hME_{cPDG}_{part}{lPDG}_any')
            any.Sumw2()
            any.Add(nonprompt)
            dME[f'hME_{cPDG}_{part}{lPDG}_any'] = any

print(dSE)
print(dME)

# project SE in the chosen pT bins
dSETmp = {}
for key in dSE:
    for ptMin, ptMax in zip(ptMins, ptMaxs):
        hh = dSE[key]
        firstBin = hh.GetYaxis().FindBin(ptMin)
        lastBin = hh.GetYaxis().FindBin(ptMax - 1.e-6)

        hProj = hh.ProjectionX(f'{key}_pt{ptMin}_{ptMax}', firstBin, lastBin)
        print(firstBin, lastBin, hProj.GetEntries())
        dSETmp[f'{key}_pt{ptMin}_{ptMax}'] = hProj
dSE = {**dSE, **dSETmp}

# project ME in the chosen pT bins
dMETmp = {}
for key in dME:
    for ptMin, ptMax in zip(ptMins, ptMaxs):
        hh = dME[key]
        firstBin = hh.GetYaxis().FindBin(ptMin)
        lastBin = hh.GetYaxis().FindBin(ptMax - 1.e-6)
        hProj = hh.ProjectionX(f'{key}_pt{ptMin}_{ptMax}', firstBin, lastBin)
        # print(firstBin, lastBin, hProj.GetEntries())

        dMETmp[f'{key}_pt{ptMin}_{ptMax}'] = hProj
dME = {**dME, **dMETmp}

cCF = TCanvas('cCF', '', 1200, 600)
cCF.Divide(2, 1)
dCF = {}


print(dME)
# Draw CF
oFile = TFile(f'{oFileNameBase}.root', 'recreate')
cCF.SaveAs(f'{oFileNameBase}_semecf.pdf[')
for lPDG in lightPDG:
    norm = GetNormRangeFromPDG(cPDG, lPDG)

    for cPDG in charmPDG:
        for part in ['', '-']:
            for ptMin, ptMax in zip(ptMins, ptMaxs):
                for charmOrigin in charmOrigins:
                    me = dME[f'hME_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}']
                    me.SetLineColor(kAzure+3)
                    me.SetLineWidth(2)
                    me.Rebin(rebin)
                    comb = 'sc' if part == '' else 'oc'
                    me.Write(f'hME_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}_pt{ptMin}_{ptMax}')
                    print(f'hME_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}')
                    print("entries me: ", me.GetEntries(), me.GetBinContent(me.GetNbinsX()+1), me.GetBinContent(me.GetNbinsX()+1)/me.GetEntries())
                    me.Scale(1./me.GetEntries())
                    me.SetTitle(';k* (GeV/c);Entries (a.u.)')
                    # me.Scale(1./me.Integral(1, me.GetNbinsX()))

                    se = dSE[f'hSE_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}']
                    se.SetLineWidth(2)
                    se.SetLineStyle(1)
                    se.SetLineColor(kRed+2)
                    se.Rebin(rebin)
                    se.Write(f'hSE_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}_pt{ptMin}_{ptMax}')
                    print("entries se: ", se.GetEntries(), se.GetBinContent(se.GetNbinsX()+1), se.GetBinContent(se.GetNbinsX()+1)/se.GetEntries())
                    se.Scale(1./se.GetEntries())
                    # se.Scale(1./se.Integral(1, se.GetNbinsX()))

                    # Draw SE/ ME distributions
                    cCF.cd(1)

                    leg = TLegend(0.4, 0.67, 0.9, 0.9)
                    leg.SetHeader(f'{ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}}')
                    leg.AddEntry(se, 'SE', 'l')
                    leg.AddEntry(me, 'ME', 'l')

                    me.GetYaxis().SetRangeUser(0, 1.5 * max([se.GetMaximum(), me.GetMaximum()]))
                    me.Draw('hist pe plc pmc')
                    se.Draw('same hist pe plc pmc')
                    leg.Draw('same')

                    # Draw correlation function
                    cCF.cd(2)
                    cf = CorrelationFunction(se=se, me=me, um='GeV2MeV', norm=norm)
                    # cf.rebin(rebin)
                    hCF = cf.get_cf()
                    if hCF is not None:
                        hCF.SetLineWidth(2)
                        hCF.GetYaxis().SetRangeUser(0.8, 1.1 * max([hCF.GetMaximum(), hCF.GetMaximum()]))
                        hCF.Write(f'hCF_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}_pt{ptMin}_{ptMax}')
                        hCF.Draw('hist pe plc pmc')
                        dCF[f'hCF_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}'] = hCF
                    else:
                        print(f"\033[33mWarning\033[0m: the CF corresponding to the key'hCF_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}' is None")
                    pairLabel = TLatex()
                    pairLabel.SetNDC()
                    pairLabel.SetTextSize(0.05)
                    pairLabel.SetTextFont(42)
                    pairLabel.DrawLatex(0.3, 0.87, TranslateToLatex(f'k{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb} {charmOrigin} D'))

                    cCF.SaveAs(f'{oFileNameBase}_semecf.pdf')
                    cCF.SaveAs(f'{oFileNameBase}_semecf_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}_pt{ptMin}_{ptMax}.png')
print('Output saved in', f'{oFileNameBase}.root')
cCF.SaveAs(f'{oFileNameBase}_semecf.pdf]')

# compare different D meson pt bins:
cCFvsPt = TCanvas("cCFvsPt", "", 600, 600)
cCFvsPt.SaveAs(f'{oFileNameBase}_ptD.pdf[')
for lPDG in lightPDG:
    for cPDG in charmPDG:
        for part in ['', '-']:
            for charmOrigin in charmOrigins:
                leg = TLegend(0.3, 0.7, 0.9, 0.9)
                comb = 'sc' if part == '' else 'oc'
                for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
                    cf = dCF[f'hCF_{cPDG}_{part}{lPDG}_{charmOrigin}_pt{ptMin}_{ptMax}']
                    if cf == None:
                        print("CF is None")
                        continue
                    if iPt == 0:
                        cf.Draw('plc pmc')
                    else:
                        cf.Draw('same plc pmc')
                    leg.AddEntry(cf, f'{ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/c')
                    leg.SetHeader(TranslateToLatex(f'k{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb} {charmOrigin} D'), 'C')
                    leg.Draw()
                cCFvsPt.SaveAs(f'{oFileNameBase}_ptD.pdf')
                cCFvsPt.SaveAs(f'{oFileNameBase}_ptD_{dPDG2Label[cPDG]}{dPDG2Label[lPDG]}_{comb}_{charmOrigin}.png')
cCFvsPt.SaveAs(f'{oFileNameBase}_ptD.pdf]')

oFile.Close()
