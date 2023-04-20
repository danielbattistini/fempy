import sys

from ROOT import gROOT, gStyle

from .correlation_function import CorrelationFunction
from .utils import *

gROOT.SetBatch(True)


# set style
def set_style(style):
    gStyle.SetPalette(55)
    gStyle.SetLineWidth(2)
    gStyle.SetLegendBorderSize(0)

    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    gStyle.SetPadRightMargin(0.035)
    gStyle.SetPadLeftMargin(0.14)
    gStyle.SetPadTopMargin(0.035)
    gStyle.SetPadBottomMargin(0.1)
    gStyle.SetOptTitle(0)

    if style == 'title':
        print('hehehe')
        gStyle.SetOptTitle(1)
        gStyle.SetPadTopMargin(0.1)

    gROOT.ForceStyle()


def error(msg):
    print(f'\033[31mError:\033[0m {msg} ---> Exit!')
    sys.exit()


def warning(msg):
    print(f'\033[33mWarning:\033[0m {msg}')


def info(msg):
    print(f'\033[34mInfo:\033[0m {msg}')
