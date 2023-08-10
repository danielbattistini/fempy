import sys

from ROOT import gROOT, gStyle

from .correlation_function import CorrelationFunction
from .utils import *

print("Welcome to fempy!")

gROOT.SetBatch(True)

import logging

class CustomFormatter(logging.Formatter):
    grey = "\x1b[90;20m"
    blue = "\x1b[34;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    purple = "\x1b[35;1m"
    reset = "\x1b[0m"
    # style = "%(module)s::%(funcName)s %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    style = "%(module)s::%(funcName)s %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + style + reset,
        logging.INFO: blue + style + reset,
        logging.WARNING: yellow + style + reset,
        logging.ERROR: red + style + reset,
        logging.CRITICAL: purple + style + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


# create logger with 'spam_application'
logger = logging.getLogger("fempy")
logger.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)

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
