'''
 --> DstarPi
python3 PlotContourYuki.py -b --pair DstarPi /home/daniel/an/DstarPi/20_luuksel/SimScanDebug_nopc_kStarBW50MeV_bs10000syst_stat_source-centr_trials.root /home/daniel/an/DstarPi/20_luuksel/SimScanDebug_nopc_kStarBW50MeV_bs10000syst_tot_source-all_trials.root /home/daniel/an/DstarPi/20_luuksel/GenCFDebug_nopc_kStarBW50MeV_bs10000syst.root /home/daniel/an/DstarPi/20_luuksel/PlotSimScanDebug_nopc_kStarBW50MeV_bs10000syst.root

 -->Pi

python3 PlotContourYuki.py -b --pair DPi /home/daniel/paper/CharmPaper/figures/final_D_20231221/SimScan_5000syst_stat_source-centr_trials.root /home/daniel/paper/CharmPaper/figures/final_D_20231221/SimScan_5000syst_tot_source-all_trials.root X /home/daniel/paper/CharmPaper/figures/final_D_20231221/PlotSimScanDebug_bs5000syst.root


'''


import argparse
import os.path as path
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from ROOT import gROOT, TFile, RDataFrame, TCanvas, TLegend, TColor, TLatex, TGraph, TGraphErrors, TLine, kGray, EColor, TGraphAsymmErrors, gStyle, gROOT
gROOT.SetBatch(True)

from fempy.utils import style


def rgb(color):
    return tuple(int(color.lstrip('#')[i:i+2], 16)/255 for i in (0, 2, 4))


def LoadPot2ScatLenConversionTable():
    dfScatLenVsPotAtr = pd.read_csv('theory/pot2a0/f0_Dpi_atr.dat', header=None, delim_whitespace=True)
    # reverse the order of the dataframe in order to avoid root getting confused with the direction of connecting lines
    dfScatLenVsPotAtr = dfScatLenVsPotAtr.iloc[::-1]

    dfScatLenVsPotRep = pd.read_csv('theory/pot2a0/f0_Dpi_rep.dat', header=None, delim_whitespace=True)
    dfScatLenVsPot = pd.concat([dfScatLenVsPotAtr, dfScatLenVsPotRep])

    gPot2ScatLen = TGraph(len(dfScatLenVsPot), dfScatLenVsPot.to_numpy()[:, 0], dfScatLenVsPot.to_numpy()[:, 1])

    return gPot2ScatLen


def plot_point_cov(points, nstd=2, ax=None, enlarge=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, enlarge, **kwargs)


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, enlarge=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    if enlarge is not None:
        shiftx, shifty = enlarge
        print(f'enlarging ellipse: {shiftx:.3f}, {shifty:.3f}')
        print(f'enlarging ellipse: {width:.3f}, {height:.3f}')

        width = np.sqrt(width**2 + 4 * shiftx**2)
        height = np.sqrt(height**2 + 4 * shifty**2)
        print(f'new: {width:.3f}, {height:.3f}')
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


def EllipseToGraph(ellipse, nPoints=500):
    gEllipse = TGraph(1)
    center = ellipse.center
    width = ellipse.width/2
    height = ellipse.height/2
    angle = ellipse.angle

    for iPoint in range(nPoints):
        x = width*np.cos(iPoint/nPoints * 2 * np.pi)
        y = height*np.sin(iPoint/nPoints * 2 * np.pi)

        xrot = + x * np.cos(angle*np.pi / 180) - y * np.sin(angle*np.pi / 180)
        yrot = + x * np.sin(angle*np.pi / 180) + y * np.cos(angle*np.pi / 180)
        
        gEllipse.SetPoint(iPoint, xrot + center[0], yrot + center[1])
    return gEllipse

albumen = '#ffdd88'
yoke = '#ff8866'

# Plot range
xlim = [-0.18, 0.38]
ylim = [-0.34, 1.099]

gROOT.SetBatch(False)

style.SetStyle(
    rightMargin=0.01,
    leftMargin=0.13,
    topMargin=0.01,
    bottomMargin=0.13)
gStyle.SetTitleOffset(1.23, "y")

parser = argparse.ArgumentParser()
parser.add_argument('inFileStat')
parser.add_argument('inFileTot')
parser.add_argument('inFileCF')
parser.add_argument('oFile')
parser.add_argument('--pair', choices=('DstarPi', 'DPi'), required=True)
parser.add_argument('-b', default=False, action='store_true')
args = parser.parse_args()

gROOT.SetBatch(args.b)

inFileTot = TFile(args.inFileTot)
tResultsTot = inFileTot.Get("tResults")
pointsTot = np.array(list(RDataFrame(tResultsTot).AsNumpy(['vq', 'vd']).values()))
dfPot = pd.DataFrame(pointsTot.T, columns=['vq', 'vd'])

# Convert potential strenght to scattering length
gPot2ScatLen = LoadPot2ScatLenConversionTable()
dfPot['a0q'] = dfPot.apply(lambda row : gPot2ScatLen.Eval(row['vq']), axis=1)
dfPot['a0d'] = dfPot.apply(lambda row : gPot2ScatLen.Eval(row['vd']), axis=1)
pointsTot = dfPot[['a0q' , 'a0d']].T.to_numpy()
_, _, xTot, yTot = dfPot.T.to_numpy()

a0qTot = xTot.mean()
a0dTot = yTot.mean()
a0qUncTot = xTot.std()
a0dUncTot = yTot.std()

# Plot a transparent 1 and 2 standard deviation covariance ellipse
ellip2Tot = plot_point_cov(pointsTot.T, nstd=2, color=albumen, label='95% CL')
ellip1Tot = plot_point_cov(pointsTot.T, nstd=1, color=yoke, label='68% CL')

cCountour = TCanvas('cCountour', '', 600, 600)
pairLatex = 'D#pi' if args.pair == 'DPi' else 'D*#pi'
titleFrame = f';#it{{a}}_{{0}}^{{{pairLatex}}} (#it{{I}} = 3/2) (fm);#it{{a}}_{{0}}^{{{pairLatex}}} (#it{{I}} = 1/2) (fm)'
hFrame = cCountour.DrawFrame(xlim[0], ylim[0], xlim[1], ylim[1], titleFrame)
hFrame.GetXaxis().CenterTitle(True)
hFrame.GetYaxis().CenterTitle(True)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)
hFrame.GetXaxis().SetNdivisions(505)
hFrame.GetYaxis().SetNdivisions(505)

tlALICE = TLatex()
tlALICE.SetTextFont(42)
tlALICE.SetTextSize(0.045)
tlALICE.DrawLatexNDC(0.23, 0.89, 'ALICE pp #sqrt{#it{s}} = 13 TeV')

tl = TLatex()
tl.SetTextFont(42)
tl.SetTextSize(0.039)
tl.DrawLatexNDC(0.23, 0.89-0.05, 'High-mult. (0#font[122]{-}0.17% INEL > 0)')

gCountourCL95 = EllipseToGraph(ellip2Tot)
gCountourCL95.SetName('gCountourCL95')
cl95Color = TColor.GetFreeColorIndex()
cl95ColorRGB = rgb(albumen)
cl95ColorRoot = TColor(cl95Color, cl95ColorRGB[0], cl95ColorRGB[1], cl95ColorRGB[2], '')
gCountourCL95.SetMarkerColor(cl95Color)
gCountourCL95.SetLineColor(cl95Color)
gCountourCL95.SetFillColor(cl95Color)
gCountourCL95.Draw('same fce3')

gCountourCL68 = EllipseToGraph(ellip1Tot)
gCountourCL68.SetName('gCountourCL68')
cl68Color = TColor.GetFreeColorIndex()
cl68RGBColor = rgb(yoke)
cl68RootColor = TColor(cl68Color, cl68RGBColor[0], cl68RGBColor[1], cl68RGBColor[2], '')
gCountourCL68.SetMarkerColor(cl68Color)
gCountourCL68.SetLineColor(cl68Color)
gCountourCL68.SetFillColorAlpha(cl68Color, 1)
gCountourCL68.Draw('same fce3')

gScattParTot = TGraphErrors(1)
gScattParTot.SetName('gScatPar')
gScattParTot.SetPoint(0, a0qTot, a0dTot)
gScattParTot.SetPointError(0, a0qUncTot, a0dUncTot)
gScattParTot.SetMarkerStyle(20)
gScattParTot.SetMarkerSize(1)
gScattParTot.SetMarkerColor(1)
gScattParTot.SetLineWidth(2)
gScattParTot.SetLineColor(1)
gScattParTot.SetMarkerColor(1)
gScattParTot.Draw('same pe')

leg = TLegend(0.48, 0.64 if args.pair == 'DstarPi' else 0.4, 0.95, 0.78)
leg.AddEntry(gScattParTot, 'Data', 'lpe')
leg.AddEntry(gCountourCL95, '95% CL', 'f')
leg.AddEntry(gCountourCL68, '68% CL', 'f')

if args.pair == 'DPi':
    lliu = TGraphAsymmErrors(1)
    lliu.SetName('gLLiu')
    lliu.SetPoint(0, -0.10, 0.37)
    lliu.SetPointError(0, 0.001, 0.001, 0.01, 0.01)
    lliu.SetMarkerColor(EColor.kRed - 7)
    lliu.SetMarkerStyle(20)
    lliu.SetMarkerSize(1.5)
    lliu.SetLineColor(EColor.kRed - 7)
    lliu.Draw('same pe')

    xyguo = TGraphErrors(1)
    xyguo.SetName('gXYGuo')
    xyguo.SetPoint(0, -0.11, 0.33)
    xyguo.SetPointError(0, 0, 0)
    xyguo.SetMarkerColor(EColor.kTeal - 9)
    xyguo.SetMarkerStyle(29)
    xyguo.SetMarkerSize(2)
    xyguo.SetLineColor(EColor.kTeal - 9)
    xyguo.Draw('same pe')

    zhguo1 = TGraphAsymmErrors(1)
    zhguo1.SetName('gZHGuo1')
    zhguo1.SetPoint(0, -0.101, 0.31)
    zhguo1.SetPointError(0, 0.005, 0.003, 0.01, 0.01)
    zhguo1.SetMarkerColor(EColor.kYellow - 6)
    zhguo1.SetMarkerStyle(34)
    zhguo1.SetMarkerSize(1.5)
    zhguo1.SetLineColor(EColor.kYellow - 6)
    zhguo1.Draw('same pe')

    zhguo2 = TGraphAsymmErrors(1)
    zhguo2.SetName('gZHGuo2')
    zhguo2.SetPoint(0, -0.099, 0.34)
    zhguo2.SetPointError(0, 0.003, 0.004, 0.00, 0.03)
    zhguo2.SetMarkerColor(EColor.kMagenta - 8)
    zhguo2.SetMarkerStyle(43)
    zhguo2.SetMarkerSize(2)
    zhguo2.SetLineColor(EColor.kMagenta - 8)
    zhguo2.Draw('same pe')

    blhuang = TGraphErrors(1)
    blhuang.SetName('gBLHuang')
    blhuang.SetPoint(0, -0.06, 0.61)
    blhuang.SetPointError(0, 0.02, 0.11)
    blhuang.SetMarkerColor(EColor.kOrange + 6)
    blhuang.SetMarkerStyle(33)
    blhuang.SetMarkerSize(1.8)
    blhuang.SetLineColor(EColor.kOrange + 6)
    blhuang.Draw('same pe')

    leg.AddEntry(xyguo, "X. Y. Guo #it{et al.}", "p")
    leg.AddEntry(zhguo1, "Z. H. Guo (Fit-1B) #it{et al.}", "p")
    leg.AddEntry(zhguo2, "Z. H. Guo (Fit-2B) #it{et al.}", "p")
    leg.AddEntry(blhuang, "B. L. Huang #it{et al.}", "p")
    leg.AddEntry(lliu, "L. Liu #it{et al.}", "p")

# same scatt par for DstarPi and DPi
gTorres = TGraphAsymmErrors(1)
gTorres.SetName('gTorres')
gTorres.SetPoint(0, -0.101, 0.423)
gTorres.SetMarkerColor(EColor.kAzure - 9)
gTorres.SetMarkerStyle(21)
gTorres.SetMarkerSize(1.5)
gTorres.SetLineColor(EColor.kAzure - 9)
gTorres.Draw('same pe')

leg.AddEntry(gTorres, "J.M. Torres-Rincon #it{et al.}", "p")
leg.Draw()

lv = TLine(0, ylim[0], 0, ylim[1]*0.7)
lv.SetLineStyle(9)
lv.SetLineColor(kGray+2)
lv.Draw()

lh = TLine(xlim[0], 0, xlim[1], 0)
lh.SetLineStyle(9)
lh.SetLineColor(kGray+2)
lh.Draw()

for ext in ['png', 'pdf', 'eps']:
    cCountour.SaveAs(f'{path.splitext(args.oFile)[0]}.{ext}')

oFile = TFile(f'{path.splitext(args.oFile)[0]}_contour.root', 'recreate')
cCountour.Write()
gCountourCL95.Write()
gCountourCL68.Write()
gScattParTot.Write()
gTorres.Write()
if (args.pair == 'DPi'):
    lliu.Write()
    xyguo.Write()
    zhguo1.Write()
    zhguo2.Write()
    blhuang.Write()

print(f'output saved in {path.splitext(args.oFile)[0]}_contour.root')

# det scat par unc
inFileStat = TFile(args.inFileStat)
tResultsStat = inFileStat.Get("tResults")
pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['vq', 'vd']).values()))
dfPotStat = pd.DataFrame(pointsStat.T, columns=['vq', 'vd'])

# Convert potential strenght to scattering length
dfPotStat['a0q'] = dfPotStat.apply(lambda row : gPot2ScatLen.Eval(row['vq']), axis=1)
dfPotStat['a0d'] = dfPotStat.apply(lambda row : gPot2ScatLen.Eval(row['vd']), axis=1)
pointsStat = dfPotStat[['a0q' , 'a0d']].T.to_numpy()
_, _, xStat, yStat = dfPotStat.T.to_numpy()

a0qStat = xStat.mean()
a0dStat = yStat.mean()
a0qUncStat = xStat.std()
a0dUncStat = yStat.std()

a0qShift = a0qStat - a0qTot
a0dShift = a0dStat - a0dTot

print()
print(f'a0q shift = {a0qShift:.8f} fm')
print(f'a0d shift = {a0dShift:.8f} fm')

a0qUncSyst = (a0qUncTot ** 2 + a0qUncStat ** 2) ** 0.5
a0dUncSyst = (a0dUncTot ** 2 + a0dUncStat ** 2) ** 0.5
print()
print(f'a0qStat = {a0qStat:.8f} +/- {a0qUncStat:.8f} (stat) +/- {a0qUncSyst:.8f} (tot ⊖ stat) fm')
print(f'a0dStat = {a0dStat:.8f} +/- {a0dUncStat:.8f} (stat) +/- {a0dUncSyst:.8f} (tot ⊖ stat) fm')

a0qUncSystShift = (a0qUncTot ** 2 + a0qUncStat ** 2 + a0qShift ** 2) ** 0.5
a0dUncSystShift = (a0dUncTot ** 2 + a0dUncStat ** 2 + a0dShift ** 2) ** 0.5

print()
print(f'a0qTot = {a0qTot:.8f} +/- {a0qUncStat:.8f} (stat) +/- {a0qUncSystShift:.8f} (tot ⊖ stat ⊕ shift) fm')
print(f'a0dTot = {a0dTot:.8f} +/- {a0dUncStat:.8f} (stat) +/- {a0dUncSystShift:.8f} (tot ⊖ stat ⊕ shift) fm')

################################################################################
# Plot range
xlim = [0, 299.99]
ylim = [0.7, 1.849]

for comb in ['SC', 'OC']:
    if args.pair == 'DstarPi':
        inFileCF = TFile(args.inFileCF)
    else:
        if comb == 'SC':
            inFileCF = TFile('/home/daniel/paper/CharmPaper/figures/final_D_20231221/PipDp_FINAL.root')
        else:
            inFileCF = TFile('/home/daniel/paper/CharmPaper/figures/final_D_20231221/PipDm_FINAL.root')

    graphsTot = inFileTot.Get(f'gCFGen{comb}')

    # graphs = GetGraphsInDir(inFileTot.Get('stat/iters'))
    # graphsSCStat = [g for g in graphs if 'SC' in g.GetName()]
    # graphsOCStat = [g for g in graphs if 'OC' in g.GetName()]

    # dgCF = {
    #     'sc': GraphAverage(graphsSCStat, graphsSCTot, shift=True),
    #     'oc': GraphAverage(graphsOCStat, graphsOCTot, shift=True),
    # }

    # # build LL curves
    # rdfTot = RDataFrame(inFileTot.Get('tot/tResults'))
    # paramsLLTot = np.array(list(rdfTot.AsNumpy(['a0sin', 'a0tri', 'r1', 'r2', 'w1']).values())).T

    # rdfStat = RDataFrame(inFileTot.Get('stat/tResults'))

    # a0sinUncSyst = np.sqrt(np.std(paramsLLTot.T[0])**2 - a0sinUncStat**2)

    # # print(f'a0(3/2) = {a0sinStat:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    # # print(f'a0(1/2) = {a0triStat:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')

    # oFile = TFile(f'{args.oFile}.root', 'create')

    # cCF = TCanvas('cCF', '', 1000, 500)
    # cCF.Divide(2, 1)
    # for iPad, comb in enumerate(['sc', 'oc']):
    #     oFile.mkdir(comb)
    #     oFile.cd(comb)

    #     pad = cCF.cd(iPad + 1)
    #     pad.SetRightMargin(0.03)
    #     pad.SetTopMargin(0.03)
    #     pad.DrawFrame(0, 0.7, 299.9999, 1.8999, TranslateToLatex(';__kStarMeV__;__C__'))


    #     # Draw syst
    #     dgCF[comb]['syst'].SetFillColorAlpha(EColor.kGray + 2, 0.65)
    #     dgCF[comb]['syst'].Draw('same e2')
    #     dgCF[comb]['syst'].Write('gCFGenSyst')
    #     # Draw stat
    #     dgCF[comb]['stat'].SetMarkerSize(1)
    #     dgCF[comb]['stat'].SetMarkerStyle(24)
    #     dgCF[comb]['stat'].Draw('same pez')
    #     dgCF[comb]['stat'].Write('gCFGenStat')

    #     # todo: uniform the calculation of the brackets with utils::ComputeBinBrackets
    #     # brackets
    #     gBrackets = TGraphErrors(1)
    #     for iPoint in range(6):
    #         gBrackets.SetPoint(iPoint, 50 * iPoint + 25, dgCF[comb]['stat'].GetPointY(iPoint))
    #         gBrackets.SetPointError(iPoint, 25, 0)

    #     SetHistStyle(gBrackets)
    #     gBrackets.SetLineWidth(1)
    #     gBrackets.SetLineColorAlpha(EColor.kBlack, 0.9)
    #     gBrackets.DrawClone("same []")
    #     gBrackets.Write('gBrackets')

    #     tl = TLatex()
    #     tl.SetTextSize(0.04)
    #     tl.SetNDC(True)
    #     if iPad == 0:
    #         tl.DrawLatex(0.22, 0.88, 'ALICE pp #sqrt{#it{s}} = 13 TeV')
    #         tl.DrawLatex(0.22, 0.82, 'High-mult. (0 #minus 0.17% INEL > 0)')
    #         tl.DrawLatex(0.22, 0.76, TranslateToLatex(f'kDstarPi_{comb}'))


    #     else:
    #         tl.DrawLatex(0.20, 0.88, TranslateToLatex(f'kDstarPi_{comb}'))
    #         a03_2 = f'a_{{0}}(I=3/2) = {a0sinStat:.2f} #pm {a0sinUncStat:.2f} (stat) #pm {a0sinUncSyst:.2f} (syst) fm'
    #         tl.DrawLatex(0.20, 0.25, a03_2)
    #         a01_2 = f'a_{{0}}(I=1/2) = {a0triStat:.2f} #pm {a0triUncStat:.2f} (stat) #pm {a0triUncSyst:.2f} (syst) fm'
    #         tl.DrawLatex(0.20, 0.25 - 0.05, a01_2)

    # for ext in ['pdf', 'png', 'eps']:
    #     cCF.SaveAs(f'{args.oFile}.{ext}')
    # oFile.cd('/')
    # cCF.Write()
















    if args.pair == 'DstarPi':
        gCFGen = inFileCF.Get(f'{comb.lower()}/gCFGenStat')
        gCFGenSystAvg = inFileCF.Get(f'{comb.lower()}/gCFGenSystAvg')
        gBrackets = inFileCF.Get(f'{comb.lower()}/gBrackets')
    else:
        gCFGen = inFileCF.Get(f'genCF_stat')
        gCFGenSystAvg = inFileCF.Get('genCF_syst')
        gBrackets = inFileCF.Get('genCF_boundaries')
    gCFFit = inFileTot.Get(f'hBestFitBand{comb}')

    # pointsTot = np.array(list(RDataFrame(tResultsTot).AsNumpy(['vq', 'vd']).values()))
    # dfPot = pd.DataFrame(pointsTot.T, columns=['vq', 'vd'])

    # # Convert potential strenght to scattering length
    # gPot2ScatLen = LoadPot2ScatLenConversionTable()
    # dfPot['a0q'] = dfPot.apply(lambda row : gPot2ScatLen.Eval(row['vq']), axis=1)
    # dfPot['a0d'] = dfPot.apply(lambda row : gPot2ScatLen.Eval(row['vd']), axis=1)
    # pointsTot = dfPot[['a0q' , 'a0d']].T.to_numpy()
    # _, _, xTot, yTot = dfPot.T.to_numpy()

    # a0qTot = xTot.mean()
    # a0dTot = yTot.mean()
    # a0qUncTot = xTot.std()
    # a0dUncTot = yTot.std()

    # # Plot a transparent 1 and 2 standard deviation covariance ellipse
    # ellip2Tot = plot_point_cov(pointsTot.T, nstd=2, color=albumen, label='95% CL')
    # ellip1Tot = plot_point_cov(pointsTot.T, nstd=1, color=yoke, label='68% CL')

    cCFGen = TCanvas(f'cCFGen{comb}', '', 600, 600)
    label = 'D'
    if args.pair == 'DstarPi':
        label += '*'
    if comb == 'SC':
        label += '^{+}#pi^{+}'
    else:
        label += '^{+}#pi^{#minus}'


    title = f';#it{{k}}* (MeV/#it{{c}});#it{{C}}_{{{label}}}(#it{{k}}*)'
    hFrame = cCountour.DrawFrame(xlim[0], ylim[0], xlim[1], ylim[1], title)
    hFrame.GetXaxis().CenterTitle(True)
    hFrame.GetYaxis().CenterTitle(True)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetLabelSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.05)
    hFrame.GetXaxis().SetNdivisions(510)
    hFrame.GetYaxis().SetNdivisions(510)

    tlALICE = TLatex()
    tlALICE.SetTextFont(42)
    tlALICE.SetTextSize(0.045)
    tlALICE.DrawLatexNDC(0.23, 0.89, 'ALICE pp #sqrt{#it{s}} = 13 TeV')

    tl = TLatex()
    tl.SetTextFont(42)
    tl.SetTextSize(0.039)
    tl.DrawLatexNDC(0.23, 0.89-0.05, 'High-mult. (0#font[122]{-}0.17% INEL > 0)')

    # gCountourCL95 = EllipseToGraph(ellip2Tot)
    # gCountourCL95.SetName('gCountourCL95')
    # cl95Color = TColor.GetFreeColorIndex()
    # cl95ColorRGB = rgb(albumen)
    # cl95ColorRoot = TColor(cl95Color, cl95ColorRGB[0], cl95ColorRGB[1], cl95ColorRGB[2], '')
    # gCountourCL95.SetMarkerColor(cl95Color)
    # gCountourCL95.SetLineColor(cl95Color)
    # gCountourCL95.SetFillColor(cl95Color)
    # gCountourCL95.Draw('same fce3')

    # gCountourCL68 = EllipseToGraph(ellip1Tot)
    # gCountourCL68.SetName('gCountourCL68')
    # cl68Color = TColor.GetFreeColorIndex()
    # cl68RGBColor = rgb(yoke)
    # cl68RootColor = TColor(cl68Color, cl68RGBColor[0], cl68RGBColor[1], cl68RGBColor[2], '')
    # gCountourCL68.SetMarkerColor(cl68Color)
    # gCountourCL68.SetLineColor(cl68Color)
    # gCountourCL68.SetFillColorAlpha(cl68Color, 1)
    # gCountourCL68.Draw('same fce3')

    # gScattParTot = TGraphErrors(1)
    # gScattParTot.SetName('gScatPar')
    # gScattParTot.SetPoint(0, a0qTot, a0dTot)
    # gScattParTot.SetPointError(0, a0qUncTot, a0dUncTot)
    # gScattParTot.SetMarkerStyle(20)
    # gScattParTot.SetMarkerSize(1)
    # gScattParTot.SetMarkerColor(1)
    # gScattParTot.SetLineWidth(2)
    # gScattParTot.SetLineColor(1)
    # gScattParTot.SetMarkerColor(1)
    # gScattParTot.Draw('same pe')

    # leg = TLegend(0.48, 0.64 if args.pair == 'DstarPi' else 0.4, 0.95, 0.78)
    # leg.AddEntry(gScattParTot, 'Data', 'lpe')
    # leg.AddEntry(gCountourCL95, '95% CL', 'f')
    # leg.AddEntry(gCountourCL68, '68% CL', 'f')

    # # same scatt par for DstarPi and DPi
    # gTorres = TGraphAsymmErrors(1)
    # gTorres.SetName('gTorres')
    # gTorres.SetPoint(0, -0.101, 0.423)
    # gTorres.SetMarkerColor(EColor.kAzure - 9)
    # gTorres.SetMarkerStyle(21)
    # gTorres.SetMarkerSize(1.5)
    # gTorres.SetLineColor(EColor.kAzure - 9)
    # gTorres.Draw('same pe')

    # leg.AddEntry(gTorres, "J.M. Torres-Rincon #it{et al.}", "p")
    # leg.Draw()

    # lv = TLine(0, ylim[0], 0, ylim[1]*0.7)
    # lv.SetLineStyle(9)
    # lv.SetLineColor(kGray+2)
    # lv.Draw()

    lh = TLine(0, 1, 300, 1)
    lh.SetLineStyle(9)
    lh.SetLineColor(kGray+2)
    lh.Draw()

    gCFFit.SetMarkerStyle(20)
    gCFFit.SetMarkerSize(0)
    gCFFit.SetFillStyle(1001)
    gCFFit.SetMarkerColorAlpha(EColor.kBlue + 1, 0.7)
    gCFFit.SetFillColorAlpha(EColor.kBlue + 1, 0.7)
    gCFFit.SetLineColorAlpha(0, 0)
    gCFFit.Draw("same e3")



    # Load Coulomb
    # for source in ['centr']:
    gCoulombCentr = inFileTot.Get(f'gCoulomb{comb}_source-centr')
    gCoulombUpper = inFileTot.Get(f'gCoulomb{comb}_source-upper')
    gCoulombLower = inFileTot.Get(f'gCoulomb{comb}_source-lower')

    gCoulomb = TGraphErrors(1)

    gCoulomb.SetName('gCOulomb')
    for iPoint in range(gCoulombCentr.GetN()):
        kStar = gCoulombCentr.GetPointX(iPoint)
        C = gCoulombCentr.GetPointY(iPoint)
        CUnc = 0.5 * abs(gCoulombUpper.GetPointY(iPoint) - gCoulombLower.GetPointY(iPoint))
        gCoulomb.SetPoint(iPoint, kStar, C)
        gCoulomb.SetPointError(iPoint, 0, CUnc)
    gCoulomb.SetFillColor(kGray + 2)
    gCoulomb.SetLineColorAlpha(0, 0)
    gCoulomb.Draw('same e3')

    leg = TLegend(0.2, 0.65, 0.8, 0.8)
    leg.SetTextSizePixels(12)
    leg.SetTextSize(0.04)
    leg.SetLineColorAlpha(0, 0)
    leg.SetFillStyle(0)
    if args.pair == 'DstarPi' and comb == 'SC':
        title = 'D*^{+}#pi^{+} #oplus D*^{#minus}#pi^{#minus}'
    elif args.pair == 'DstarPi' and comb == 'OC':
        title = 'D*^{+}#pi^{#minus} #oplus D*^{#minus}#pi^{+}'

    elif args.pair == 'DPi' and comb == 'SC':
        title = 'D^{+}#pi^{+} #oplus D^{#minus}#pi^{#minus}'
    elif args.pair == 'DPi' and comb == 'OC':
        title = 'D^{+}#pi^{#minus} #oplus D^{#minus}#pi^{+}'
        
    leg.AddEntry(gCFGen, title, 'pel')
    leg.AddEntry(gCoulomb, 'Coulomb', 'f')
    leg.AddEntry(gCFFit, 'Gaussian potential', 'f')
    leg.Draw('same')




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
        
    # SetHistStyle(gBrackets)
    for iPoint in range(gBrackets.GetN()):
        if type(gCFGenSystAvg) == TGraphErrors:
            gCFGenSystAvg.SetPointError(iPoint, 0.5 * gCFGen.GetErrorX(iPoint), gCFGenSystAvg.GetErrorY(iPoint))
        elif type(gCFGenSystAvg) == TGraphAsymmErrors:
            gCFGenSystAvg.SetPointError(iPoint, 0.5 * gCFGen.GetErrorX(iPoint), 0.5 * gCFGen.GetErrorX(iPoint), gCFGenSystAvg.GetErrorY(iPoint), gCFGenSystAvg.GetErrorY(iPoint))

    gCFGenSystAvg.SetLineColorAlpha(0, 0)
    gCFGenSystAvg.SetFillColorAlpha(kGray+1, 0.65)
    gCFGenSystAvg.SetFillStyle(1001)
    gCFGenSystAvg.Draw('same ze2')

    # gCFGen.SetMarkerStyle(24)
    gCFGen.SetMarkerStyle(20)
    gCFGen.SetMarkerSize(0.5)
    gCFGen.Draw('same zpe')

    for iPoint in range(gBrackets.GetN()):
        gBrackets.SetPoint(iPoint, gBrackets.GetPointX(iPoint), gCFGen.GetPointY(iPoint))
    
    gBrackets.SetLineWidth(1)
    gBrackets.SetLineColor(EColor.kBlack)
    gBrackets.Draw("same []")
    

    if comb == 'OC':
        pair = 'D#pi' if args.pair == 'DPi' else 'D*#pi'
        a0qStat = round(a0qStat, 2)
        if a0qStat <= 0:
            if a0qStat > 0.01:
                a0qStat = f'{abs(a0qStat):.2f}'
            else:
                a0qStat = f'{a0qStat:.2f}'.replace('-', '#minus')
        else:
            a0qStat = f'{a0qStat:.2f}'

        a0dStat = round(a0dStat, 2)
        if a0dStat <= 0:
            if a0dStat > -0.01:
                a0dStat = f'{abs(a0dStat):.2f}'
            else:
                a0dStat = f'{a0dStat:.2f}'.replace('-', '#minus')
        else:
            a0dStat = f'{a0dStat:.2f}'

        a03_2 = f'a^{{{pair}}}_{{0}}(I=3/2) = {a0qStat} #pm {a0qUncStat:.2f} (stat.) #pm {a0qUncSystShift:.2f} (syst.) fm'
        tl.DrawLatexNDC(0.20, 0.25, a03_2)
        a01_2 = f'a^{{{pair}}}_{{0}}(I=1/2) = {a0dStat} #pm {a0dUncStat:.2f} (stat.) #pm {a0dUncSystShift:.2f} (syst.) fm'
        tl.DrawLatexNDC(0.20, 0.25 - 0.05, a01_2)

    for ext in ['png', 'pdf', 'eps']:
        cCFGen.SaveAs(f'{path.splitext(args.oFile)[0]}_fit{comb}.{ext}')

    oFile = TFile(f'{path.splitext(args.oFile)[0]}_fit{comb}.root', 'recreate')

oFile.Close()
