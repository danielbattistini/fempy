import sys
import argparse
import os.path as path

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from ROOT import gROOT, TFile, RDataFrame, TCanvas, TLegend, TColor, TLatex, TGraph, TGraphErrors, TLine, kGray, EColor, TGraphAsymmErrors, gStyle
gROOT.SetBatch(True)

def rgb(color):
    return tuple(int(color.lstrip('#')[i:i+2], 16)/255 for i in (0, 2, 4))


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
    gStyle.SetPadColor(10)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetHistLineWidth(1)
    gStyle.SetFuncWidth(2)
    gStyle.SetLineWidth(2)
    gStyle.SetLabelSize(0.045, "xyz")
    # gStyle.SetLabelOffset(0.01, "y")
    # gStyle.SetLabelOffset(0.01, "x")
    gStyle.SetTitleSize(0.05, "xyz")
    gStyle.SetTitleOffset(1.15, "y")
    # gStyle.SetTitleOffset(1.2, "x")
    gStyle.SetTextSizePixels(26)
    gStyle.SetTextFont(42)
    gStyle.SetLegendFont(42)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetErrorX(0.005)

# DstarPi
# python3 PlotContour.py /home/daniel/an/DstarPi/20_luuksel/SimFit_nopc_kStarBW50MeV_bs10002syst.root /home/daniel/an/DstarPi/20_luuksel/SimFitContour_DstarPi_nopc_kStarBW50MeV_bs10002syst_stat_and_shift_final.png  --central statbs --shift --pair DstarPi
# python3 PlotContour.py /home/daniel/an/DstarPi/20_luuksel/SimFit_nopc_kStarBW50MeV_bs10002syst_source_and_fit_range_var_as_syst_Dstarmass_scalableLL.root /home/daniel/an/DstarPi/20_luuksel/SimFitContour_nopc_kStarBW50MeV_bs10002syst_source_and_fit_range_var_as_syst_Dstarmass_scalableLL.root --shift --pair DstarPi --central statbs
# python3 PlotContour.py /home/daniel/an/DstarPi/20_luuksel/SimFit_nopc_kStarBW50MeV_bs10003syst_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700_chi2ndflt5_originalinputfile_indepRandGen.root SimFitContour_nopc_kStarBW50MeV_bs10003syst_uncThermalFist-beauty-DstarPurity_fixQSRedMasSwapp_combfitLL_scaledLL_fit700_chi2ndflt5_originalinputfile_indepRandGen.root --shift --pair DstarPi --central statbs

# DPi
# python3 PlotContour.py /home/daniel/an/DPi/simfit/SimFit.root /home/daniel/an/DPi/simfit/SimFitContourPlot_DPi_check2.png  --central statbs --shift --pair DPi --emma
# python3 PlotContour.py /home/daniel/an/DPi/simfit/SimFit2.root /home/daniel/an/DPi/simfit/SimFitContourPlot_DPi_scalableLL.png  --central statbs --shift --pair DPi
# python3 PlotContour.py /home/daniel/an/DPi/simfit/SimFit2.root /home/daniel/an/DPi/simfit/SimFitContourPlot_DPi_scalableLL_yuki  --central statbs --shift --pair DPi --yuki ~/paper/CharmPaper/figures/final_D/simscan/SimScan_SimScan_fast_check3_trials.root

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('oFile')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--shift', action='store_true', default=False)
    parser.add_argument('--emma', action='store_true', default=False)
    parser.add_argument('--central', choices=('statbs', 'totbs'))
    parser.add_argument('--pair', choices=('DstarPi', 'DPi'), required=True)
    parser.add_argument('--yuki', default=None, help='file with the trials of the CF a la yuki')
    args = parser.parse_args()

    Setstyle()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    inFile = TFile(args.inFile)

    tResultsTot = inFile.Get("tot/tResults")
    if args.emma:
        pointsTot = np.array(list(RDataFrame(tResultsTot).AsNumpy(['a0tri', 'a0sin']).values()))
    else:
        pointsTot = np.array(list(RDataFrame(tResultsTot).AsNumpy(['a0sin', 'a0tri']).values()))
    xTot, yTot = pointsTot
    a0sinTot = xTot.mean()
    a0triTot = yTot.mean()
    a0sinUncTot = xTot.std()
    a0triUncTot = yTot.std()
    albumen = '#ffdd88'
    yoke = '#ff8866'

    # Load stat bootstrap
    tResultsStat = inFile.Get("stat/tResults")
    if args.emma:
        pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0tri', 'a0sin']).values()))
    else:
        pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))
    xStat, yStat = pointsStat
    a0sinStat = xStat.mean()
    a0triStat = yStat.mean()
    a0sinUncStat = xStat.std()
    a0triUncStat = yStat.std()

    # Load results a la yuki
    if args.yuki:
        inFileYuki = TFile(args.yuki)
        tYuki = inFileYuki.Get('tResults')
        pointsYuki = np.array(list(RDataFrame(tYuki).AsNumpy(['a0sin', 'a0tri']).values()))
    
    # Make figure
    xlim = [-0.18, 0.38]
    ylim = [-0.34, 1.099]

    # Plot a transparent 3 standard deviation covariance ellipse
    ellip2Tot = plot_point_cov(pointsTot.T, nstd=2, color=albumen, label='95% CL')
    ellip1Tot = plot_point_cov(pointsTot.T, nstd=1, color=yoke, label='68% CL')
    
    ellip2TotShifted = plot_point_cov(pointsTot.T, nstd=2, color=albumen, label='95% CL')
    ellip2TotShifted.center = (a0sinStat, a0triStat)
    ellip1TotShifted = plot_point_cov(pointsTot.T, nstd=1, color=yoke, label='68% CL')
    ellip1TotShifted.center = (a0sinStat, a0triStat)

    absShift = (abs(a0sinTot - a0sinStat), abs(a0triTot - a0triStat))
    ellip2TotShiftedEnlarged = plot_point_cov(pointsTot.T, nstd=2, color=albumen, enlarge=absShift, label='95% CL')
    ellip2TotShiftedEnlarged.center = (a0sinStat, a0triStat)
    ellip1TotShiftedEnlarged = plot_point_cov(pointsTot.T, nstd=1, color=yoke, enlarge=absShift, label='68% CL')
    ellip1TotShiftedEnlarged.center = (a0sinStat, a0triStat)

    ellip2Stat = plot_point_cov(pointsStat.T, nstd=2, color=albumen, label='95% CL')
    ellip1Stat = plot_point_cov(pointsStat.T, nstd=1, color=yoke, label='68% CL')

    a0sinUncSyst = (a0sinUncTot**2 - a0sinUncStat**2)**0.5
    a0triUncSyst = (a0triUncTot**2 - a0triUncStat**2)**0.5

    if args.debug:
        pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))

    cCountour = TCanvas('cCountour', '', 600, 600)
    cCountour.SetRightMargin(0.01)
    cCountour.SetTopMargin(0.01)
    cCountour.SetLeftMargin(0.13)
    cCountour.SetBottomMargin(0.13)
    pairLatex = 'D#pi' if args.pair == 'DPi' else 'D*#pi'
    hFrame = cCountour.DrawFrame(xlim[0], ylim[0], xlim[1], ylim[1], f';#it{{a}}_{{0}}^{{{pairLatex}}} (#it{{I}} = 3/2) (fm);#it{{a}}_{{0}}^{{{pairLatex}}} (#it{{I}} = 1/2) (fm)')
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

    tl = TLatex()
    tl.SetTextFont(42)
    tl.SetTextSize(0.039)

    gScattParTot = TGraphErrors(1)
    gScattParTot.SetPoint(0, a0sinTot, a0triTot)
    gScattParTot.SetPointError(0, a0sinUncTot, a0triUncTot)
    gScattParTot.SetMarkerStyle(20)
    gScattParTot.SetMarkerSize(1)
    gScattParTot.SetMarkerColor(1)
    gScattParTot.SetLineWidth(2)

    gScattParStat = TGraphErrors(1)
    gScattParStat.SetPoint(0, a0sinStat, a0triStat)
    gScattParStat.SetPointError(0, a0sinUncTot, a0triUncTot)
    gScattParStat.SetMarkerStyle(20)
    gScattParStat.SetMarkerSize(1)
    gScattParStat.SetMarkerColor(1)
    gScattParStat.SetLineWidth(2)

    gScattParStatCentrTotPlusShift = TGraphErrors(1)
    gScattParStatCentrTotPlusShift.SetPoint(0, a0sinStat, a0triStat)
    a0sinUncTotPlusShift = np.sqrt(a0sinUncTot**2 + absShift[0]**2)
    a0triUncTotPlusShift = np.sqrt(a0triUncTot**2 + absShift[1]**2)
    a0sinUncSystPlusShift = np.sqrt(a0sinUncTotPlusShift**2 - a0sinUncStat**2)
    a0triUncSystPlusShift = np.sqrt(a0triUncTotPlusShift**2 - a0triUncStat**2)
    gScattParStatCentrTotPlusShift.SetPointError(0, a0sinUncTotPlusShift, a0triUncTotPlusShift)

    print('\ncentr total')
    print(f'a0(3/2) = {a0sinTot:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    print(f'a0(1/2) = {a0triTot:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')

    print('\ncentr stat')
    print(f'a0(3/2) = {a0sinStat:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    print(f'a0(1/2) = {a0triStat:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')

    print('\ncentr stat, shift and enlarge')
    print(f'a0(3/2) = {a0sinStat:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSystPlusShift:.2f} (syst) fm')
    print(f'a0(1/2) = {a0triStat:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSystPlusShift:.2f} (syst) fm')

    gCountour2Tot = EllipseToGraph(ellip2Tot)
    gCountour1Tot = EllipseToGraph(ellip1Tot)
    gCountour2Stat = EllipseToGraph(ellip2Stat)
    gCountour1Stat = EllipseToGraph(ellip1Stat)

    gCountour1TotShifted = EllipseToGraph(ellip1TotShifted)
    gCountour2TotShifted = EllipseToGraph(ellip2TotShifted)

    gCountour1TotShiftedEnlarged = EllipseToGraph(ellip1TotShiftedEnlarged)
    gCountour2TotShiftedEnlarged = EllipseToGraph(ellip2TotShiftedEnlarged)

    if args.central == 'statbs':
        if args.shift:
            gContourToFill1 = gCountour1TotShiftedEnlarged
            gContourToFill2 = gCountour2TotShiftedEnlarged
        else:
            gContourToFill1 = gCountour1Stat
            gContourToFill2 = gCountour2Stat

        gContourToDash1 = gCountour1Tot
        gContourToDash2 = gCountour2Tot
            
        dLeg = {
            'fill1': '68% CL stat bs',
            'fill2': '95% CL stat bs',
            'dash1': '68% CL tot bs',
            'dash2': '95% CL tot bs',
        }
    elif args.central == 'totbs':
        gContourToFill1 = gCountour1Tot
        gContourToFill2 = gCountour2Tot
    
        gContourToDash1 = gCountour1Stat
        gContourToDash2 = gCountour2Stat
        dLeg = {
            'fill1': '68% CL tot bs',
            'fill2': '95% CL tot bs',
            'dash1': '68% CL stat bs',
            'dash2': '95% CL stat bs',
        }

    ciAlbumen = TColor.GetFreeColorIndex()
    albumenRGB = rgb(albumen)
    albumenRoot = TColor(ciAlbumen, albumenRGB[0], albumenRGB[1], albumenRGB[2], '')
    gContourToFill2.SetFillColor(ciAlbumen)

    yokeRGB = rgb(yoke)
    ciYoke = TColor.GetFreeColorIndex()
    yokeRoot = TColor(ciYoke, yokeRGB[0], yokeRGB[1], yokeRGB[2], '')
    gContourToFill1.SetFillColorAlpha(ciYoke, 1)

    gContourToDash1.SetLineStyle(9)
    gContourToDash1.SetLineWidth(2)
    gContourToDash1.SetLineColor(4)

    gContourToDash2.SetLineStyle(9)
    gContourToDash2.SetLineWidth(2)
    gContourToDash2.SetLineColor(4)

    lv = TLine(0, ylim[0], 0, ylim[1]*0.7)

    lv.SetLineStyle(9)
    lv.SetLineColor(kGray+2)

    lh = TLine(xlim[0], 0, xlim[1], 0)
    lh.SetLineStyle(9)
    lh.SetLineColor(kGray+2)
    
    if args.pair == 'DstarPi':
        legPosX = [0.48, 0.95]
        legPosY = [0.64, 0.78]
    elif args.pair == 'DPi':
        legPosX = [0.48, 0.95]
        legPosY = [0.4, 0.78]
    leg = TLegend(legPosX[0], legPosY[0], legPosX[1], legPosY[1])
    if args.debug:
        leg = TLegend(0.34, 0.55, 0.95, 0.95)
        leg.AddEntry(gScattParTot, 'Data tot bs', 'lpe')
        if args.central == 'statbs':
            gScattParStat.SetLineColor(1)
            gScattParStat.SetMarkerColor(1)
            gScattParTot.SetLineColor(4)
            gScattParTot.SetMarkerColor(4)
        else:
            gScattParStat.SetLineColor(4)
            gScattParStat.SetMarkerColor(4)
            gScattParTot.SetLineColor(1)
            gScattParTot.SetMarkerColor(1)

        leg.AddEntry(gScattParStat, 'Data stat bs', 'lpe')
        leg.AddEntry(gScattParStatCentrTotPlusShift, 'Data stat bs, tot #oplus shift', 'lpe')

        leg.AddEntry(gContourToFill1, dLeg['fill1'], 'f')
        leg.AddEntry(gContourToFill2, dLeg['fill2'], 'f')
        leg.AddEntry(gContourToDash1, dLeg['dash1'], 'l')
        leg.AddEntry(gContourToDash2, dLeg['dash2'], 'l')
        gContourToFill2.Draw('same fce3')
        gContourToFill1.Draw('same fce3')
        gContourToDash2.Draw('same l')
        gContourToDash1.Draw('same l')

        gCountour2TotShifted.SetLineStyle(8)
        gCountour2TotShifted.SetLineColor(2)

        gCountour1TotShifted.SetLineStyle(8)
        gCountour1TotShifted.SetLineColor(2)
        gCountour1TotShifted.Draw('same l')
        gCountour2TotShifted.Draw('same l')

        leg.AddEntry(gCountour1TotShifted, '68% CL tot bs shifted', 'l')
        leg.AddEntry(gCountour2TotShifted, '95% CL tot bs shifted', 'l')

        gCountour2TotShiftedEnlarged.SetLineStyle(7)
        gCountour2TotShiftedEnlarged.SetLineColor(3)

        gCountour1TotShiftedEnlarged.SetLineStyle(7)
        gCountour1TotShiftedEnlarged.SetLineColor(3)
        gCountour1TotShiftedEnlarged.Draw('same l')
        gCountour2TotShiftedEnlarged.Draw('same l')

        leg.AddEntry(gCountour1TotShiftedEnlarged, '68% CL tot bs shifted enlarged', 'l')
        leg.AddEntry(gCountour2TotShiftedEnlarged, '95% CL tot bs shifted enlarged', 'l')
        
        gScattParStatCentrTotPlusShift.SetMarkerStyle(25)
        gScattParStatCentrTotPlusShift.SetMarkerSize(2)
        gScattParStatCentrTotPlusShift.SetMarkerColor(617)
        gScattParStatCentrTotPlusShift.SetLineColor(617)
        gScattParStatCentrTotPlusShift.SetLineWidth(2)
        gScattParStatCentrTotPlusShift.Draw('same lpe')

        gScattParTot.Draw('same pe')

    else:
        leg.AddEntry(gContourToFill1, '68% CL', 'f')
        leg.AddEntry(gContourToFill2, '95% CL', 'f')
        
        gContourToFill2.Draw('same fce3')
        gContourToFill1.Draw('same fce3')
            
        leg.AddEntry(gScattParStatCentrTotPlusShift, 'Data', 'lpe')

        gCountour1TotShiftedEnlarged.SetLineColorAlpha(0, 0)
        gCountour2TotShiftedEnlarged.SetLineColorAlpha(0, 0)
        gCountour1TotShiftedEnlarged.Draw('same l')
        gCountour2TotShiftedEnlarged.Draw('same l')

        if args.yuki:
            ellip2Yuki = plot_point_cov(pointsYuki.T, nstd=2, label='95% CL')
            ellip1Yuki = plot_point_cov(pointsYuki.T, nstd=1, label='68% CL')
            gCountour2Yuki = EllipseToGraph(ellip2Yuki)
            gCountour1Yuki = EllipseToGraph(ellip1Yuki)

            gCountour2Yuki.SetLineStyle(7)
            gCountour1Yuki.SetLineStyle(7)
            gCountour2Yuki.SetLineColor(EColor.kAzure+2)
            gCountour1Yuki.SetLineColor(EColor.kAzure+2)
            gCountour2Yuki.SetLineWidth(1)
            gCountour1Yuki.SetLineWidth(1)
            
            gCountour2Yuki.Draw('same')
            gCountour1Yuki.Draw('same')
            leg.AddEntry(gCountour1Yuki, 'Gaussian potential', 'l')

        gScattParStatCentrTotPlusShift.SetMarkerStyle(20)
        gScattParStatCentrTotPlusShift.SetMarkerSize(1)
        gScattParStatCentrTotPlusShift.SetMarkerColor(1)
        gScattParStatCentrTotPlusShift.SetLineColor(1)
        gScattParStatCentrTotPlusShift.SetLineWidth(2)
        gScattParStatCentrTotPlusShift.Draw('same lpe')

    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)

    lv.Draw()
    lh.Draw()
    leg.Draw()

    if not args.debug:
        tlALICE.DrawLatexNDC(0.23, 0.89, 'ALICE pp #sqrt{#it{s}} = 13 TeV')
        tl.DrawLatexNDC(0.23, 0.89-0.05, 'High-mult. (0#font[122]{-}0.17% INEL > 0)')

    if args.debug:
        gScattParStat.Draw('same pez')

    if args.pair == 'DPi':
        lliu = TGraphAsymmErrors(1)
        lliu.SetPoint(0, -0.10, 0.37)
        lliu.SetPointError(0, 0.001, 0.001, 0.01, 0.01)
        lliu.SetMarkerColor(EColor.kRed - 7)
        lliu.SetMarkerStyle(20)
        lliu.SetMarkerSize(1.5)
        lliu.SetLineColor(EColor.kRed - 7)

        xyguo = TGraphErrors(1)
        xyguo.SetPoint(0, -0.11, 0.33)
        xyguo.SetPointError(0, 0, 0)
        xyguo.SetMarkerColor(EColor.kTeal - 9)
        xyguo.SetMarkerStyle(29)
        xyguo.SetMarkerSize(2)
        xyguo.SetLineColor(EColor.kTeal - 9)

        zhguo1 = TGraphAsymmErrors(1)
        zhguo1.SetPoint(0, -0.101, 0.31)
        zhguo1.SetPointError(0, 0.005, 0.003, 0.01, 0.01)
        zhguo1.SetMarkerColor(EColor.kYellow - 6)
        zhguo1.SetMarkerStyle(34)
        zhguo1.SetMarkerSize(1.5)
        zhguo1.SetLineColor(EColor.kYellow - 6)

        zhguo2 = TGraphAsymmErrors(1)
        zhguo2.SetPoint(0, -0.099, 0.34)
        zhguo2.SetPointError(0, 0.003, 0.004, 0.00, 0.03)
        zhguo2.SetMarkerColor(EColor.kMagenta - 8)
        zhguo2.SetMarkerStyle(43)
        zhguo2.SetMarkerSize(2)
        zhguo2.SetLineColor(EColor.kMagenta - 8)

        blhuang = TGraphErrors(1)
        blhuang.SetPoint(0, -0.06, 0.61)
        blhuang.SetPointError(0, 0.02, 0.11)
        blhuang.SetMarkerColor(EColor.kOrange + 6)
        blhuang.SetMarkerStyle(33)
        blhuang.SetMarkerSize(1.8)
        blhuang.SetLineColor(EColor.kOrange + 6)

        lliu.Draw('same pe')
        xyguo.Draw('same pe')
        zhguo1.Draw('same pe')
        zhguo2.Draw('same pe')
        blhuang.Draw('same pe')

        leg.AddEntry(xyguo, "X. Y. Guo #it{et al.}", "p")
        leg.AddEntry(zhguo1, "Z. H. Guo (Fit-1B) #it{et al.}", "p")
        leg.AddEntry(zhguo2, "Z. H. Guo (Fit-2B) #it{et al.}", "p")
        leg.AddEntry(blhuang, "B. L. Huang #it{et al.}", "p")
        leg.AddEntry(lliu, "L. Liu #it{et al.}", "p")

    # same scatt par for DstarPi and DPi
    gTorres = TGraphAsymmErrors(1)
    gTorres.SetPoint(0, -0.101, 0.423)
    gTorres.SetMarkerColor(EColor.kAzure - 9)
    gTorres.SetMarkerStyle(21)
    gTorres.SetMarkerSize(1.5)
    gTorres.SetLineColor(EColor.kAzure - 9)
    leg.AddEntry(gTorres, "J.M. Torres-Rincon #it{et al.}", "p")
    gTorres.Draw('same pe')

    for ext in ['png', 'pdf', 'eps', 'root']:
        if args.debug:
            name = f'{path.splitext(args.oFile)[0]}_debug.{ext}'
        else:
            name = f'{path.splitext(args.oFile)[0]}.{ext}'
        cCountour.SaveAs(name)
