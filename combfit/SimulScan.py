import pandas as pd
import numpy as np
import argparse

from ROOT import TFile, TCanvas, TGraph, TF1, TLegend, SetOwnership, gStyle, TH2D, kRed, kGray, TLatex, gROOT, EColor

gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
# default = '~/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs10002syst.root'
# parser.add_argument('inFile', help='input file without extension')
parser.add_argument('oFileWoExt', help='output file without extension')
parser.add_argument('--pair', choices=('DPi', 'DstarPi'), help='DPi or DstarPi')
args = parser.parse_args()


sourceRadii = [0.97, 2.52]
weightSource = 0.66
baseDir = '../theory/cf/yuki/DPi_correlation'

gCFTheoSC = {}
gCFTheoOC = {}

fitRange = [0, 300]

for vq in range(-1000, 2000, 100):
    try:
        dfsc_small = pd.read_csv(f'{baseDir}/fast_check2/corr_Dp_pip/0.97fm/corr_0.97fm_{vq}.dat', header=None, delim_whitespace=True)
        dfsc_large = pd.read_csv(f'{baseDir}/fast_check2/corr_Dp_pip/2.52fm/corr_2.52fm_{vq}.dat', header=None, delim_whitespace=True)
    except FileNotFoundError:
        print(f'File not found for vq = {vq}')
        continue
    dfsc = dfsc_small
    dfsc[1] = weightSource * dfsc_small[1] + (1 - weightSource) * dfsc_large[1]
    gCFTheoSC[vq] = TGraph(len(dfsc), dfsc.to_numpy()[:, 0], dfsc.to_numpy()[:, 1])
    gCFTheoSC[vq].SetName(f'gCFTheoSC_Vq{vq}')

    gCFTheoOC[vq] = {}
    # Load opposite charge CFs
    for vd in range(-1000, 2000, 100):
        try:
            dfsc_small = pd.read_csv(f'{baseDir}/fast_check2/corr_Dp_pim/0.97fm/corr_0.97fm_Vq_{vq}_Vd_{vd}.dat', header=None, delim_whitespace=True)
            dfsc_large = pd.read_csv(f'{baseDir}/fast_check2/corr_Dp_pim/2.52fm/corr_2.52fm_Vq_{vq}_Vd_{vd}.dat', header=None, delim_whitespace=True)
        except FileNotFoundError:
            print(f'File not found for vq = {vq} and vd = {vd}')
            continue
        dfsc = dfsc_small
        dfsc[1] = weightSource * dfsc_small[1] + (1 - weightSource) * dfsc_large[1]
        gCFTheoOC[vq][vd] = TGraph(len(dfsc), dfsc.to_numpy()[:, 0], dfsc.to_numpy()[:, 1])
        gCFTheoOC[vq][vd].SetName(f'gCFTheoOC_Vq{vq}_Vd{vd}')

# Plot the loaded CFs
gStyle.SetPalette(55)
cCFTheoSC = TCanvas('cCFTheoSC', '', 600, 600)
frame = cCFTheoSC.DrawFrame(0, fitRange[0], fitRange[1], 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
leg = TLegend(0.15, 0.5, 0.85, 0.85)
leg.SetNColumns(2)
for vq, fCFSC in gCFTheoSC.items():
    leg.AddEntry(fCFSC, f'Vq={vq} MeV')
    fCFSC.Draw('pl same plc pmc')
leg.Draw()
cCFTheoSC.SaveAs(f'{args.oFileWoExt}_theoSC.pdf')

cCFTheoOC = TCanvas('cCFTheoOC', '', 600, 600)
cCFTheoOC.SaveAs(f'{args.oFileWoExt}_theoOC.pdf[')
for vq, _ in gCFTheoOC.items():
    frame = cCFTheoOC.DrawFrame(0, fitRange[0], fitRange[1], 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
    leg = TLegend(0.15, 0.5, 0.85, 0.85, f'Vq = {vq} MeV')
    leg.SetNColumns(2)
    for vd in gCFTheoOC[vq].keys():
        leg.AddEntry(gCFTheoOC[vq][vd], f'Vd={vd} MeV')
        gCFTheoOC[vq][vd].Draw('pl same plc pmc')
    leg.Draw()

    cCFTheoOC.SaveAs(f'{args.oFileWoExt}_theoOC.pdf')
cCFTheoOC.SaveAs(f'{args.oFileWoExt}_theoOC.pdf]')

# Load data
def LoadCFGen(pair):
    if pair == 'DPi':
        inFileSC = TFile('/home/daniel/paper/CharmPaper/figures/final_D/PipDp_FINAL.root')
        gsc = inFileSC.Get('genCF_stat')
        inFileOC = TFile('/home/daniel/paper/CharmPaper/figures/final_D/PipDm_FINAL.root')
        goc = inFileOC.Get('genCF_stat')
    elif pair == 'DstarPi':
        inFile = TFile('~/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs10002syst.root')
        gsc = inFile.Get('sc/stat/gCFGen0')
        goc = inFile.Get('oc/stat/gCFGen0')
    return(gsc, goc)

gCFGenSC, gCFGenOC = LoadCFGen(args.pair)

def ComputeChi2(data, theory, xmin, xmax):
    chi2 = 0
    for iPoint in range(data.GetN()):
        x = data.GetPointX(iPoint)
        if x < xmin:
            continue
        if x > xmax:
            break

        y = data.GetPointY(iPoint)
        yUnc = data.GetErrorY(iPoint)
        yTheo = theory.Eval(x)

        chi2 += ((y - yTheo)/yUnc) ** 2

    return chi2

cFits = TCanvas('cFits', '', 1200, 600)
cFits.Divide(2, 1)
cFits.SaveAs(f'{args.oFileWoExt}_fits.pdf[')
hhChi2 = TH2D('hhChi2', '', 30, -1000, 2000, 30, -1000, 2000)
hhChi2Ndf = TH2D('hhChi2Ndf', '', 30, -1000, 2000, 30, -1000, 2000)

if args.pair == 'DPi':
    titleSC = 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}'
    titleOC = 'D^{+}#pi^{#minus} #oplus D^{#minus}#pi^{+}'
elif args.pair == 'DstarPi':
    titleSC = 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}'
    titleOC = 'D*^{+}#pi^{#minus} #oplus D*^{#minus}#pi^{+}'

minChi2 = 1.e10
bestPot = (None, None)


def CountPoints(graph, xmin, xmax):
    count = 0

    for iPoint in range(graph.GetN()):
        x = graph.GetPointX(iPoint)
        if xmin < x and x < xmax:
            count += 1
    return count

ndf = CountPoints(gCFGenSC, fitRange[0], fitRange[1]) + CountPoints(gCFGenOC, fitRange[0], fitRange[1]) - 2

for vq, gSC in gCFTheoSC.items():
    for vd, gOC in gCFTheoOC[vq].items():
        pad = cFits.cd(1)
        frame = pad.DrawFrame(0, 0, 300, 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
        gSC.SetMarkerColor(kRed)
        gSC.SetLineColor(kRed)
        gSC.Draw('same l')
        gCFGenSC.Draw('same pe')
        gCFTheoSC[0].SetMarkerColor(kGray+2)
        gCFTheoSC[0].SetLineColor(kGray+2)
        gCFTheoSC[0].SetLineStyle(9)
        gCFTheoSC[0].Draw('same l')
        chi2SC = ComputeChi2(gCFGenSC, gSC, fitRange[0], fitRange[1])

        legSC = TLegend(0.3, 0.7, 0.88, 0.88, f'Vq = {vq} MeV')


        legSC.AddEntry(gCFGenSC, titleSC)
        legSC.AddEntry(gSC, f'theory (#chi^{{2}} = {chi2SC:.1f})')
        legSC.AddEntry(gCFTheoSC[0], f'Coulomb (#chi^{{2}} = {chi2SC:.1f})')
        legSC.Draw()

        pad = cFits.cd(2)
        frame = pad.DrawFrame(0, fitRange[0], fitRange[1], 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
        gOC.SetLineColor(kRed)
        chi2OC = ComputeChi2(gCFGenOC, gOC, fitRange[0], fitRange[1])
        gOC.Draw('same l')
        gCFGenOC.Draw('same pe')

        gCFTheoOC[0][0].SetMarkerColor(kGray+2)
        gCFTheoOC[0][0].SetLineColor(kGray+2)
        gCFTheoOC[0][0].SetLineStyle(9)
        gCFTheoOC[0][0].Draw('same l')

        legOC = TLegend(0.3, 0.7, 0.88, 0.88, f'Vq = {vq} MeV; Vd = {vd} MeV')
        legOC.AddEntry(gCFGenOC, titleOC)
        legOC.AddEntry(gOC, f'theory (#chi^{{2}} = {chi2OC:.1f})')
        legOC.AddEntry(gCFTheoOC[0][0], f'Coulomb (#chi^{{2}} = {chi2OC:.1f})')
        legOC.Draw()

        cFits.SaveAs(f'{args.oFileWoExt}_fits.pdf')

        binx = hhChi2.GetXaxis().FindBin(vq)
        biny = hhChi2.GetYaxis().FindBin(vd)
        chi2 = chi2SC + chi2OC
        hhChi2.SetBinContent(binx, biny, chi2)

        hhChi2Ndf.SetBinContent(binx, biny, chi2/ndf)

        if chi2 < minChi2:
            minChi2 = chi2
            bestPot = (vq, vd)

cFits.SaveAs(f'{args.oFileWoExt}_fits.pdf]')

cChi2 = TCanvas('cChi2', '', 600, 600)
cChi2.SetLeftMargin(0.15)
cChi2.SetRightMargin(0.18)
gStyle.SetOptStat(0)
hhChi2.SetTitle(';Vq (MeV);Vd (MeV);#chi^{2}')
hhChi2.Draw('same colz')

gBestPot = TGraph(1)
gBestPot.SetPoint(0, bestPot[0], bestPot[1])
gBestPot.SetMarkerStyle(20)
gBestPot.SetMarkerSize(1)
gBestPot.SetMarkerColor(EColor.kOrange+10)
gBestPot.SetLineColor(EColor.kOrange+10)
gBestPot.SetLineColor(EColor.kOrange+10)
gBestPot.Draw('pe same')

cChi2.SaveAs(f'{args.oFileWoExt}_contour.pdf')



cChi2Ndf = TCanvas('cChi2Ndf', '', 600, 600)
cChi2Ndf.SetLeftMargin(0.15)
cChi2Ndf.SetRightMargin(0.18)
gStyle.SetOptStat(0)
hhChi2Ndf.SetTitle(';Vq (MeV);Vd (MeV);#chi^{2}/ndf')
hhChi2Ndf.Draw('same colz')
hContour = hhChi2Ndf.Clone('hContour')
hContour.SetContour(2, np.array([(minChi2 + 2.30)/ndf, (minChi2 + 6.18) / ndf], 'd'))
hContour.SetLineColor(EColor.kRed)
hContour.Draw('same cont3')


gBestPot = TGraph(1)
print(bestPot, minChi2, minChi2/ndf)
gBestPot.SetPoint(0, bestPot[0], bestPot[1])
gBestPot.SetMarkerStyle(20)
gBestPot.SetMarkerSize(1)
gBestPot.SetMarkerColor(EColor.kOrange+10)
gBestPot.SetLineColor(EColor.kOrange+10)
gBestPot.SetLineColor(EColor.kOrange+10)
gBestPot.Draw('pe same')

cChi2Ndf.SaveAs(f'{args.oFileWoExt}_Chi2Ndf.pdf')
