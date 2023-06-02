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
    fig, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    xlim = [-0.07, 0.2999]
    ylim = [-0.2999, 0.3999]
    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.text(0.23, 0.9, 'ALICE', fontsize=22, transform=ax.transAxes)
    plt.text(0.23, 0.9-0.07, 'pp $\sqrt{s}$ = 13 TeV, High-mult. (0-0.17%)', fontsize=13, transform=ax.transAxes)
    plt.plot([0, 0], ylim, '--', color='gray')
    plt.plot(xlim, [0,0], '--', color='gray')

    # Plot a transparent 3 standard deviation covariance ellipse
    ellip2 = plot_point_cov(pointsTot.T, nstd=2, color=albumen, label='95% CL')
    ellip1 = plot_point_cov(pointsTot.T, nstd=1, color=yoke, label='68% CL')
    plt.legend(handles=[ellip1, ellip2], loc='upper left', bbox_to_anchor=(0.6,0.8), fontsize=13, frameon=False)
    plt.errorbar([a0sinTot], [a0triTot], xerr=a0sinUncTot, yerr=a0triUncTot, color='black', marker='o')


    a0sinUncSyst = (a0sinUncTot**2 - a0sinUncStat**2)**0.5
    a0triUncSyst = (a0triUncTot**2 - a0triUncStat**2)**0.5

    print(f'a0(3/2) = {a0sinTot:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    print(f'a0(1/2) = {a0triTot:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')

    plt.xlabel('$a_{0}~(I=3/2)$', fontsize=18)
    plt.ylabel('$a_{0}~(I=1/2)$', fontsize=18)
    fig.tight_layout()

    if args.debug:
        pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))
        plt.errorbar([a0sinStat], [a0triStat], xerr=a0sinUncStat, yerr=a0triUncStat, color='blue', marker='o')
        ellip2Stat = plot_point_cov(
            pointsStat.T,
            nstd=2,
            fill=False,
            edgecolor='blue',
            linestyle='--',
            label='95% CL stat')
        ellip1Stat = plot_point_cov(
            pointsStat.T,
            nstd=1,
            fill=False,
            edgecolor='blue',
            linestyle='--',
            label='68% CL stat')
        plt.legend(
            handles=[ellip1, ellip2, ellip1Stat, ellip2Stat],
            loc='upper left',
            bbox_to_anchor=(0.6,0.8),
            fontsize=13,
            frameon=False)
    
    plt.show()

    plt.savefig(f'{path.splitext(args.oFile)[0]}_py_debug{path.splitext(args.oFile)[1]}' if args.debug else f'{path.splitext(args.oFile)[0]}_py{path.splitext(args.oFile)[1]}')
 
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
    tlALICE.DrawLatexNDC(0.36, 0.89, 'ALICE')

    tl = TLatex()
    tl.SetTextFont(42)
    tl.SetTextSize(0.039)
    tl.DrawLatexNDC(0.36, 0.89-0.05, 'pp #sqrt{s} = 13 TeV, High-mult. (0-0.17%)')

    gScattPar = TGraphErrors(1)

    gScattPar.SetPoint(0, a0sinTot, a0triTot)
    gScattPar.SetPointError(0, a0sinUncTot, a0triUncTot)
    gScattPar.SetMarkerStyle(20)
    gScattPar.SetMarkerSize(1)
    gScattPar.SetMarkerColor(1)
    gScattPar.SetLineWidth(2)

    gCountour2 = EllipseToGraph(ellip2)
    ciAlbumen = TColor.GetFreeColorIndex()
    albumenRGB = rgb(albumen)
    albumenRoot = TColor(ciAlbumen, albumenRGB[0], albumenRGB[1], albumenRGB[2], '')
    gCountour2.SetFillColor(ciAlbumen)

    gCountour1 = EllipseToGraph(ellip1)
    yokeRGB = rgb(yoke)
    ciYoke = TColor.GetFreeColorIndex()
    yokeRoot = TColor(ciYoke, yokeRGB[0], yokeRGB[1], yokeRGB[2], '')
    gCountour1.SetFillColorAlpha(ciYoke, 1)

    lh = TLine(0, ylim[0], 0, ylim[1])
    lh.SetLineStyle(9)
    lh.SetLineColor(kGray+2)

    lv = TLine(xlim[0], 0, xlim[1], 0)
    lv.SetLineStyle(9)
    lv.SetLineColor(kGray+2)
    
    if args.debug:
        leg = TLegend(0.34, 0.57, 0.95, 0.78)
    else:
        leg = TLegend(0.65, 0.57, 0.95, 0.78)
    leg.AddEntry(gScattPar, 'Data', 'lpe')
    leg.AddEntry(gCountour1, '68% CL', 'f')
    leg.AddEntry(gCountour2, '95% CL', 'f')
    leg.SetBorderSize(0)
    leg.Draw()

    gCountour2.Draw('same fce3')
    gCountour1.Draw('same fce3')
    lv.Draw()
    lh.Draw()
    gScattPar.Draw('same pe')
    cCountour.SaveAs(f'{path.splitext(args.oFile)[0]}_root_debug{path.splitext(args.oFile)[1]}'if args.debug else f'{path.splitext(args.oFile)[0]}_root{path.splitext(args.oFile)[1]}')
