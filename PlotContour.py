import sys
import argparse
import os.path as path

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


from ROOT import TFile, RDataFrame


def plot_point_cov(points, nstd=2, ax=None, **kwargs):
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
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
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
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('oFile')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    
    inFile = TFile(args.inFile)

    tResults = inFile.Get("tot/tResults")
    points = np.array(list(RDataFrame(tResults).AsNumpy(['a0sin', 'a0tri']).values()))
    x, y = points
    a0sin = x.mean()
    a0tri = y.mean()
    a0sinUncTot = x.std()
    a0triUncTot = y.std()

    points = np.array(list(RDataFrame(tResults).AsNumpy(['a0sin', 'a0tri']).values()))
    x, y = points

    fig, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    xlim = [-0.07, 0.3]
    ylim = [-0.4, 0.4]
    plt.xlim(xlim)
    plt.ylim(ylim)

    step = 0.07
    plt.text(0.23, 0.9, 'ALICE', fontsize=22, transform=ax.transAxes)
    plt.text(0.23, 0.9-step, 'pp $\sqrt{s}$ = 13 TeV, High-mult. (0-0.17%)', fontsize=13, transform=ax.transAxes)
    plt.plot([0, 0], ylim, '--', color='gray')
    plt.plot(xlim, [0,0], '--', color='gray')

    # Plot a transparent 3 standard deviation covariance ellipse
    ellip3 = plot_point_cov(points.T, nstd=3, color='#ffdd88', label='95% CL')
    ellip1 = plot_point_cov(points.T, nstd=1, color='orange', label='68% CL')
    plt.legend(handles=[ellip1, ellip3], loc='upper left', bbox_to_anchor=(0.6,0.8), fontsize=13, frameon=False)
    plt.errorbar([a0sin], [a0tri], xerr=a0sinUncTot, yerr=a0triUncTot, color='black', marker='o',  markeredgecolor="black", markerfacecolor="black")


    # stat unc
    tResultsStat = inFile.Get("stat/tResults")
    points = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))
    x, y = points
    a0sinStat = x.mean()
    a0triStat = y.mean()
    a0sinUncStat = x.std()
    a0triUncStat = y.std()

    a0sinUncSyst = (a0sinUncTot**2 - a0sinUncStat**2)**0.5
    a0triUncSyst = (a0triUncTot**2 - a0triUncStat**2)**0.5

    print(f'a0(3/2) = {a0sin:.2f} +/- {a0sinUncStat:.2f} (stat) +/- {a0sinUncSyst:.2f} (syst) fm')
    print(f'a0(1/2) = {a0tri:.2f} +/- {a0triUncStat:.2f} (stat) +/- {a0triUncSyst:.2f} (syst) fm')
    
    plt.xlabel('$a_0~(I=3/2)$', fontsize=18)
    plt.ylabel('$a_0~(I=1/2)$', fontsize=18)
    fig.tight_layout()

    if args.debug:

        pointsStat = np.array(list(RDataFrame(tResultsStat).AsNumpy(['a0sin', 'a0tri']).values()))
        plt.errorbar([a0sinStat], [a0triStat], xerr=a0sinUncStat, yerr=a0triUncStat, color='blue', facecolors=None, marker='o',  markeredgecolor="blue", markerfacecolor="blue")
        ellip3Stat = plot_point_cov(pointsStat.T, nstd=3, fill=False, edgecolor='blue',linestyle='--', label='95% CL stat')
        ellip1Stat = plot_point_cov(pointsStat.T, nstd=1, fill=False, edgecolor='blue',linestyle='--', label='68% CL stat')
        plt.legend(handles=[ellip1, ellip3, ellip1Stat, ellip3Stat], loc='upper left', bbox_to_anchor=(0.6,0.8), fontsize=13, frameon=False)
    

    plt.show()
    plt.savefig(f'{path.splitext(args.oFile)[0]}_debug{path.splitext(args.oFile)[1]}'if args.debug else args.oFile)
    
    
    
    
    
