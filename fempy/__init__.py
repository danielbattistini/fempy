import sys
import logging

from ROOT import gROOT, gStyle

from . import utils
from . import sim

# Create logging formatter for beautiful log messages
class CustomFormatter(logging.Formatter):
    '''
    DEBUG: for very detailed output
    INFO: for tracking the progress of the program
    WARNING: for something unexpected that doesn't cause problems
    ERROR: for something unexpected that can cause problems
    CRITICAL: for a failure that cause the program to terminate
    '''
    grey = "\x1b[90;20m"
    blue = "\x1b[34;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    purple = "\x1b[35;1m"
    reset = "\x1b[0m"
    style = "%(module)s::%(funcName)s %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    criticalStyle = "%(module)s::%(funcName)s %(levelname)s - %(message)s ---> EXIT! (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + style + reset,
        logging.INFO: blue + style + reset,
        logging.WARNING: yellow + style + reset,
        logging.ERROR: red + style + reset,
        logging.CRITICAL: purple + criticalStyle + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

# Create handler to exit the program for CRITICAL errors
class ExitOnExceptionHandler(logging.StreamHandler):
    '''Exit on CRITICAL'''

    def emit(self, record):
        super().emit(record)
        if record.levelno is logging.CRITICAL:
            raise SystemExit(-1)


print("Welcome to fempy!")
gROOT.SetBatch(True)

# create logger
logger = logging.getLogger("fempy")
logger.setLevel(logging.INFO)
ch = ExitOnExceptionHandler()
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
