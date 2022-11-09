from ROOT import gROOT, gStyle

from .correlation_function import CorrelationFunction

gROOT.SetBatch(True)

# set style
# gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(55)
gStyle.SetLineWidth(2)
gStyle.SetLegendBorderSize(0)

gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetPadBottomMargin(0.1)
