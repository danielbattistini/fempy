import sys
import argparse
import os.path as path

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from ROOT import TFile, RDataFrame, TCanvas, TLegend, TColor, TLatex, TGraph, TGraphErrors, TLine, kGray

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
        width = np.sqrt(width**2 + 4 * shiftx**2)
        height = np.sqrt(height**2 + 4 * shifty**2)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('oFile')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--shift', action='store_true', default=False)
    parser.add_argument('--central', choices=('statbs', 'totbs'))
    args = parser.parse_args()

    inFile = TFile(args.inFile)

    tResultsTot = inFile.Get("tot/tResults")
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
    pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))
    xStat, yStat = pointsStat
    a0sinStat = xStat.mean()
    a0triStat = yStat.mean()
    a0sinUncStat = xStat.std()
    a0triUncStat = yStat.std()

    # Make figure
    xlim = [-0.07, 0.2999]
    ylim = [-0.2999, 0.3999]

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
    cCountour.SetRightMargin(0.03)
    cCountour.SetTopMargin(0.03)
    cCountour.SetLeftMargin(0.18)
    cCountour.SetBottomMargin(0.12)
    hFrame = cCountour.DrawFrame(xlim[0], ylim[0], xlim[1], ylim[1], ';a_{0} (#it{I} = 3/2);a_{0} (#it{I} = 1/2)')
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetLabelSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.05)
    hFrame.GetXaxis().SetNdivisions(504)
    hFrame.GetYaxis().SetNdivisions(504)

    tlALICE = TLatex()
    tlALICE.SetTextFont(42)
    tlALICE.SetTextSize(0.06)
    if not args.debug:
        tlALICE.DrawLatexNDC(0.36, 0.89, 'ALICE')

    tl = TLatex()
    tl.SetTextFont(42)
    tl.SetTextSize(0.039)
    if not args.debug:
        tl.DrawLatexNDC(0.36, 0.89-0.05, 'pp #sqrt{s} = 13 TeV, High-mult. (0-0.17%)')

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

    lh = TLine(0, ylim[0], 0, ylim[1])
    lh.SetLineStyle(9)
    lh.SetLineColor(kGray+2)

    lv = TLine(xlim[0], 0, xlim[1], 0)
    lv.SetLineStyle(9)
    lv.SetLineColor(kGray+2)
    
    leg = TLegend(0.56, 0.57, 0.95, 0.78)
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

        gScattParStatCentrTotPlusShift.SetMarkerStyle(20)
        gScattParStatCentrTotPlusShift.SetMarkerSize(1)
        gScattParStatCentrTotPlusShift.SetMarkerColor(1)
        gScattParStatCentrTotPlusShift.SetLineColor(1)
        gScattParStatCentrTotPlusShift.SetLineWidth(2)
        gScattParStatCentrTotPlusShift.Draw('same lpe')

    leg.SetBorderSize(0)

    lv.Draw()
    lh.Draw()
    leg.Draw()

    if args.debug:
        gScattParStat.Draw('same pe')

    for ext in ['png', 'pdf', 'eps', 'root']:
        if args.debug:
            name = f'{path.splitext(args.oFile)[0]}_root_debug.{ext}'
        else:
            name = f'{path.splitext(args.oFile)[0]}_root.{ext}'
        cCountour.SaveAs(name)
