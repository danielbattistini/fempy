import sys
import pandas as pd
import numpy as np
import argparse

from ROOT import TFile, TCanvas, TGraph, TLegend, TNtuple, gStyle, TH2D, gROOT, EColor, TProfile, TGraphErrors

gROOT.SetBatch(True)


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


def CountPoints(graph, xmin, xmax):
    count = 0

    for iPoint in range(graph.GetN()):
        x = graph.GetPointX(iPoint)
        if xmin < x and x < xmax:
            count += 1
    return count


parser = argparse.ArgumentParser()
# default = '~/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs10002syst.root'
# parser.add_argument('inFile', help='input file without extension')
parser.add_argument('oFileWoExt', help='output file without extension')
parser.add_argument('--unc', choices=('stat', 'tot'), default='stat', help='Stat or Syst run')
parser.add_argument('--pair', choices=('DPi', 'DstarPi'), help='DPi or DstarPi')
parser.add_argument('--source', choices=('centr', 'upper', 'lower'), help='Source variation', default='centr')
args = parser.parse_args()

if args.source == 'centr':
    sourceRadii = [0.97, 2.52]
    weightSource = 0.66
elif args.source == 'lower':
    sourceRadii = [0.89, 2.32]
    weightSource = 0.69
elif args.source == 'upper':
    sourceRadii = [1.06, 2.88]
    weightSource = 0.64
    
baseDir = '../theory/cf/yuki'

gCFTheoSC = {}
gCFTheoOC = {}

fitRanges = [[0, 250], [0, 200], [0, 300]]

for vq in range(-3000, 2100, 100):
    try:
        print(f'{baseDir}/Dp_Pip/{sourceRadii[0]}fm/corr_{sourceRadii[0]}fm_{vq}.dat')
        dfsc_small = pd.read_csv(f'{baseDir}/Dp_Pip/{sourceRadii[0]}fm/corr_{sourceRadii[0]}fm_{vq}.dat', header=None, delim_whitespace=True)
        dfsc_large = pd.read_csv(f'{baseDir}/Dp_Pip/{sourceRadii[1]}fm/corr_{sourceRadii[1]}fm_{vq}.dat', header=None, delim_whitespace=True)
    except FileNotFoundError:
        print(f'File not found for vq = {vq}')
        continue
    dfsc = dfsc_small
    dfsc[1] = weightSource * dfsc_small[1] + (1 - weightSource) * dfsc_large[1]
    gCFTheoSC[vq] = TGraph(len(dfsc), dfsc.to_numpy()[:, 0], dfsc.to_numpy()[:, 1])
    gCFTheoSC[vq].SetName(f'gCFTheoSC_Vq{vq}')

    gCFTheoOC[vq] = {}
    # Load opposite charge CFs
    for vd in range(-2000, 5100, 100):
        try:
            dfsc_small = pd.read_csv(f'{baseDir}/Dp_Pim/{sourceRadii[0]}fm/corr_{sourceRadii[0]}fm_Vq_{vq}_Vd_{vd}.dat', header=None, delim_whitespace=True)
            dfsc_large = pd.read_csv(f'{baseDir}/Dp_Pim/{sourceRadii[1]}fm/corr_{sourceRadii[1]}fm_Vq_{vq}_Vd_{vd}.dat', header=None, delim_whitespace=True)
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
frame = cCFTheoSC.DrawFrame(0, 0, 300, 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
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
    frame = cCFTheoOC.DrawFrame(0, 0, 300, 3, ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
    leg = TLegend(0.15, 0.5, 0.85, 0.85, f'Vq = {vq} MeV')
    leg.SetNColumns(2)
    for vd in gCFTheoOC[vq].keys():
        leg.AddEntry(gCFTheoOC[vq][vd], f'Vd={vd} MeV')
        gCFTheoOC[vq][vd].Draw('pl same plc pmc')
    leg.Draw()

    cCFTheoOC.SaveAs(f'{args.oFileWoExt}_theoOC.pdf')
cCFTheoOC.SaveAs(f'{args.oFileWoExt}_theoOC.pdf]')

# Load data

cFits = TCanvas('cFits', '', 1200, 600)
cFits.Divide(2, 1)
cFits.SaveAs(f'{args.oFileWoExt}_fits.pdf[')

nBinsVq = len(gCFTheoSC)
minVq = min(gCFTheoSC.keys())
maxVq = max(gCFTheoSC.keys())
print(gCFTheoSC.keys())
stepVq = round((maxVq - minVq)/(nBinsVq -1))

nBinsVd = len(gCFTheoOC[0])
minVd = min(gCFTheoOC[0].keys())
maxVd = max(gCFTheoOC[0].keys())
stepVd = round((maxVd - minVd)/(nBinsVd -1))

print(nBinsVq, minVq, maxVq, stepVq)
print(nBinsVd, minVd, maxVd, stepVd)


hhChi2 = TH2D('hhChi2', '', nBinsVq, minVq-stepVq/2, maxVq+stepVq/2, nBinsVd, minVd-stepVd/2, maxVd+stepVd/2)
hhChi2Ndf = TH2D('hhChi2Ndf', '', nBinsVq, minVq-stepVq/2, maxVq+stepVq/2, nBinsVd, minVd-stepVd/2, maxVd+stepVd/2)


if args.pair == 'DPi':
    titleSC = 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}'
    titleOC = 'D^{+}#pi^{#minus} #oplus D^{#minus}#pi^{+}'
elif args.pair == 'DstarPi':
    titleSC = 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}'
    titleOC = 'D*^{+}#pi^{#minus} #oplus D*^{#minus}#pi^{+}'

# Load conversion table between potentials and scatt len
dfScatLenVsPotAtr = pd.read_csv('../theory/pot2a0/f0_Dpi_atr.dat', header=None, delim_whitespace=True)
# reverse the order of the dataframe in order to avoid root getting confused with  direction of connecting lines
dfScatLenVsPotAtr = dfScatLenVsPotAtr.iloc[::-1]

dfScatLenVsPotRep = pd.read_csv('../theory/pot2a0/f0_Dpi_rep.dat', header=None, delim_whitespace=True)
dfScatLenVsPot = pd.concat([dfScatLenVsPotAtr, dfScatLenVsPotRep])

gPot2ScatLen = TGraph(len(dfScatLenVsPot), dfScatLenVsPot.to_numpy()[:, 0], dfScatLenVsPot.to_numpy()[:, 1])

gStyle.SetOptStat(0)

tResults = TNtuple("tResults", "", "vq:vd:chi2ndf:a0sin:a0tri")

cChi2Ndf = TCanvas('cChi2Ndf', '', 1800, 600)
cChi2Ndf.Divide(3, 1)
pad = cChi2Ndf.cd(1)
pad.SetLeftMargin(0.15)
pad.SetRightMargin(0.02)
pad = cChi2Ndf.cd(2)
pad.SetLeftMargin(0.15)
pad.SetRightMargin(0.02)
pad = cChi2Ndf.cd(3)
pad.SetLeftMargin(0.15)
pad.SetRightMargin(0.21)
cChi2Ndf.SaveAs(f'{args.oFileWoExt}_Chi2Ndf.pdf[')

if args.pair == 'DPi':
    # inFileSC = TFile("/home/daniel/paper/CharmPaper/figures/final_D/simscan/DPiPlusOutput_tot.root")
    # inFileOC = TFile("/home/daniel/paper/CharmPaper/figures/final_D/simscan/DPiMinusOutput_tot.root")
    inFileSC = TFile("/home/daniel/paper/CharmPaper/figures/final_D_20231221/DPiPlusOutput_tot.root")
    inFileOC = TFile("/home/daniel/paper/CharmPaper/figures/final_D_20231221/DPiMinusOutput_tot.root")
    gCFGenSC = inFileSC.Get("CF_corr_MC_truth_NBL_0")
    gCFGenOC = inFileOC.Get("CF_corr_MC_truth_NBL_0")
elif args.pair == 'DstarPi':
    # inFile = TFile('~/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs10002syst.root')
    # inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs1000_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700_chi2ndflt5_originalinputfile_indepRandGen_normrange300-600.root')
    inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs50_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700_chi2ndflt5_originalinputfile_indepRandGen_normrange1500-2000_bkgfitrange300-1000.root')
    inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs5000syst_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_noLLfit_originalinputfile_indepRandGen_normrange1500-2000_bkgfitrange300-1000.root')
    # inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/GenCFCorrTest_nopc_kStarBW50MeV_bs100syst_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700_chi2ndflt5_originalinputfile_indepRandGen.root')
    # inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/GenCFDebug_nopc_kStarBW50MeV_bs10000syst.root')

    gCFGenSC = inFile.Get(f'sc/{args.unc}/gCFGen0')
    gCFGenOC = inFile.Get(f'oc/{args.unc}/gCFGen0')

gCoulombSC = gCFTheoSC[0]
gCoulombOC = gCFTheoOC[0][0]
xvals = gCoulombSC.GetX()

# Assume that all bins have a width of 1 MeV
edges = np.array(xvals, 'd') + 0.5
edges = np.insert(edges, 0, 0.5, axis=0)
hBestFitBandSC = TProfile("hBestFitBandSC", "", gCFTheoSC[0].GetN(), edges, 's')
hBestFitBandOC = TProfile("hBestFitBandOC", "", gCFTheoSC[0].GetN(), edges, 's')

oFile = TFile(f"{args.oFileWoExt}_trials.root", 'create')
if oFile.IsZombie(): # File already exists
    sys.exit()


for iIter in range(100000):
    print(f"Iter: {iIter}")

    fitRange = fitRanges[np.random.randint(len(fitRanges)) if args.unc == 'tot' else 0]
    ndf = CountPoints(gCFGenSC, fitRange[0], fitRange[1]) + CountPoints(gCFGenOC, fitRange[0], fitRange[1]) - 2

    if args.pair == 'DPi':
        gCFGenSC = inFileSC.Get(f"CF_corr_MC_truth_NBL_{iIter}")
        gCFGenOC = inFileOC.Get(f"CF_corr_MC_truth_NBL_{iIter}")
    elif args.pair == 'DstarPi':
        gCFGenSC = inFile.Get(f'sc/{args.unc}/gCFGen{iIter}')
        gCFGenOC = inFile.Get(f'oc/{args.unc}/gCFGen{iIter}')

    if gCFGenSC == None or gCFGenOC == None:
        break

    # Start minimization
    minChi2 = 1.e10
    bestPot = (None, None)

    for vq, gSC in gCFTheoSC.items():
        chi2SC = ComputeChi2(gCFGenSC, gSC, fitRange[0], fitRange[1])
        for vd, gOC in gCFTheoOC[vq].items():
            chi2OC = ComputeChi2(gCFGenOC, gOC, fitRange[0], fitRange[1])
            chi2 = chi2SC + chi2OC

            binx = hhChi2.GetXaxis().FindBin(vq)
            biny = hhChi2.GetYaxis().FindBin(vd)
            hhChi2Ndf.SetBinContent(binx, biny, chi2/ndf)

            if chi2 < minChi2:
                minChi2 = chi2
                bestPot = (vq, vd)
    for iPoint in range(gCFTheoSC[bestPot[0]].GetN()):
        xx = gCFTheoSC[bestPot[0]].GetPointX(iPoint)
        yy = gCFTheoSC[bestPot[0]].GetPointY(iPoint)
        hBestFitBandSC.Fill(xx, yy)
        
    for iPoint in range(gCFTheoOC[bestPot[0]][bestPot[1]].GetN()):
        xx = gCFTheoOC[bestPot[0]][bestPot[1]].GetPointX(iPoint)
        yy = gCFTheoOC[bestPot[0]][bestPot[1]].GetPointY(iPoint)
        hBestFitBandOC.Fill(xx, yy)

    # print(bestPot[0], bestPot[1], gCFTheoOC[bestPot[0]][bestPot[1]].GetPointY(30), gCFTheoOC[0][0].GetPointY(30))
    tResults.Fill(bestPot[0], bestPot[1], minChi2 / ndf, gPot2ScatLen.Eval(bestPot[0]), gPot2ScatLen.Eval(bestPot[1]))

    if iIter < 50:
        # Draw best-fit CF (same-charge)
        pad = cChi2Ndf.cd(1)
        frame = pad.DrawFrame(0, 0.7, 300, 1.3, 'same-charge;#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
        gCFTheoSC[0].SetLineColor(EColor.kGray+2)
        gCFTheoSC[0].SetLineStyle(7)
        gCFTheoSC[0].Draw('same')
        gCFTheoSC[bestPot[0]].SetLineColor(EColor.kAzure + 2)
        gCFTheoSC[bestPot[0]].Draw('same')
        gCFGenSC.Draw('same')
        legScanSC = TLegend(0.4, 0.7, 0.9, 0.9)
        legScanSC.AddEntry(gCFGenSC, 'Data')
        legScanSC.AddEntry(gCFTheoSC[0], f'Coulomb', 'l')
        legScanSC.AddEntry(gCFTheoSC[bestPot[0]], f'Best fit: V_{{q}} = {bestPot[0]} MeV', 'l')
        legScanSC.Draw('same')

        # Draw best-fit CF (opposite-charge)
        pad = cChi2Ndf.cd(2)
        frame = pad.DrawFrame(0, 0.7, 300, 1.3, 'opposite-charge;#it{k}* (MeV/#it{c});#it{C}(#it{k}*)')
        gCFTheoOC[0][0].SetLineColor(EColor.kGray+2)
        gCFTheoOC[0][0].SetLineStyle(7)
        gCFTheoOC[0][0].Draw('same')
        gCFTheoOC[bestPot[0]][bestPot[1]].SetLineColor(EColor.kAzure + 2)
        gCFTheoOC[bestPot[0]][bestPot[1]].Draw('same')
        gCFGenOC.Draw('same')
        legScanOC = TLegend(0.4, 0.7, 0.9, 0.9)
        legScanOC.AddEntry(gCFGenOC, 'Data')
        legScanOC.AddEntry(gCFTheoOC[0][0], f'Coulomb', 'l')
        legScanOC.AddEntry(gCFTheoOC[bestPot[0]][bestPot[1]], f'Best fit: V_{{q}} = {bestPot[0]} MeV, V_{{d}} = {bestPot[1]} MeV', 'l')
        legScanOC.Draw('same')

        # Draw Chi2/ndf plot
        pad = cChi2Ndf.cd(3)
        hhChi2Ndf.SetTitle(f'iIter = {iIter};Vq (MeV);Vd (MeV);#chi^{{2}}/ndf')
        hhChi2Ndf.SetNdivisions(510, 'xy')
        hhChi2Ndf.Draw('colz')
        hContour = hhChi2Ndf.Clone('hContour')
        hContour.SetContour(2, np.array([(minChi2 + 2.30) / ndf, (minChi2 + 6.18) / ndf], 'd'))
        hContour.SetLineColor(EColor.kRed)
        hContour.Draw('same cont3')

        gBestPot = TGraph(1)
        # print(bestPot, minChi2, minChi2/ndf)
        gBestPot.SetPoint(0, bestPot[0], bestPot[1])
        gBestPot.SetMarkerStyle(20)
        gBestPot.SetMarkerSize(1)
        gBestPot.SetMarkerColor(EColor.kOrange+10)
        gBestPot.SetLineColor(EColor.kOrange+10)
        gBestPot.SetLineColor(EColor.kOrange+10)
        gBestPot.Draw('p same')

        cChi2Ndf.Modified()
        cChi2Ndf.Update()
        cChi2Ndf.SaveAs(f'{args.oFileWoExt}_Chi2Ndf.pdf')

cChi2Ndf.SaveAs(f'{args.oFileWoExt}_Chi2Ndf.pdf]')

gCoulombSC.SetName(f'gCoulombSC_source-{args.source}')
gCoulombSC.Write()

gCoulombOC.SetName(f'gCoulombOC_source-{args.source}')
gCoulombOC.Write()

hBestFitBandSC.Write()
hBestFitBandOC.Write()

def SubtractInQuadrature(graph1, graph2, name='', shift=True):
    if graph1.GetN() != graph2.GetN():
        print("diff num of points")
        sys.exit()

    gAns = TGraphErrors(1)
    gAns.SetName(name)
    for iPoint in range(graph1.GetN()):
        gAns.SetPoint(iPoint, graph1.GetPointX(iPoint), graph1.GetPointY(iPoint))

        e1 = graph1.GetErrorY(iPoint)
        e2 = graph2.GetErrorY(iPoint)
        systSq = e1 ** 2 - e2 ** 2
        if shift:
            y1 = graph1.GetErrorY(iPoint)
            y2 = graph2.GetErrorY(iPoint)
            systSq += (y1 - y2) ** 2

        if systSq < 0:
            print("warning: Syst unc < 0. Assigning 0")
            systSq = 0

        gAns.SetPointError(iPoint, graph1.GetErrorX(iPoint), systSq ** 0.5)

    return gAns

if args.pair == 'DstarPi':
    gCFGenSCstat = inFile.Get(f'sc/gCFGenStat')
    gCFGenSCstat.SetName('gCFGenSCstat')
    gCFGenSCsyst = inFile.Get(f'sc/gCFGenSyst')
    gCFGenSCsyst.SetName('gCFGenSCsyst')

    gCFGenOCstat = inFile.Get(f'oc/gCFGenStat')
    gCFGenOCstat.SetName('gCFGenOCstat')
    gCFGenOCsyst = inFile.Get(f'oc/gCFGenSyst')
    gCFGenOCsyst.SetName('gCFGenOCsyst')

else:
    inFileGenCF_SC = TFile('/home/daniel/paper/CharmPaper/figures/final_D_20231221/PipDp_FINAL.root')
    gCFGenSCstat = inFileGenCF_SC.Get('genCF_stat')
    gCFGenSCstat.SetName('gCFGenSCstat')
    gCFGenSCsyst = inFileGenCF_SC.Get('genCF_syst')
    gCFGenSCsyst.SetName('gCFGenSCsyst')

    inFileGenCF_OC = TFile('/home/daniel/paper/CharmPaper/figures/final_D_20231221/PipDm_FINAL.root')
    gCFGenOCstat = inFileGenCF_OC.Get('genCF_stat')
    gCFGenOCstat.SetName('gCFGenOCstat')
    gCFGenOCsyst = inFileGenCF_OC.Get('genCF_syst')
    gCFGenOCsyst.SetName('gCFGenOCsyst')

oFile.cd()
gCFGenSCstat.Write()
gCFGenSCsyst.Write()
gCFGenOCstat.Write()
gCFGenOCsyst.Write()

tResults.Write()
oFile.Close()
print(f'Output saved in {args.oFileWoExt}_trials.root')
