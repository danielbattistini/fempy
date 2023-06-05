'''
Script to create a tree dataset from binned dataset
'''
import numpy as np
import argparse
from sys import exit

from utils.io import GetObjectFromFile
from utils.handle import GetFemtoDreamPairId
from utils.correlation_handler import CorrelationHandler

from ROOT import TFile, TTree, gRandom, TH1D, TKDE, RDataFrame, TF1, SetOwnership, gROOT, TCanvas

# def get(se, me, x):
#     return se.Eval(x) / me.Eval(x)

num, den = TF1(), TF1()

def extract_raw_cf(inFileName, oFileName):
    """
    Extract the raw correlation function from the raw CF, sidebands, minijets etc. CFs.

    Parameters
    ----------
    inFileName : str
        path of input datafile

    Returns
    -------
    int
        Description of return value

    """

    pairs = ['sc']
    regions = ['sgn', 'sbl', 'sbr']
    kdeEstimates = {}
    CorrelationFunctions = {}

    # oFile = TFile(oFileName, 'recreate')
    # for pairId in pairs:
    #     for regionId in regions:

    def compute_ratio_root(x):
        return num.Eval(x[0]) / den.Eval(x[0])

    in_file = TFile(inFileName)
    num = in_file.Get('SE_sc_sgn_centr')
    den = in_file.Get('ME_sc_sgn_centr')

    getCorrRoot = lambda x, par: compute_ratio_root(x)
    # SetOwnership(se, False)
    # se.SetName('SE_sc_sgn_centr')
    # se.AddToGlobalList()

    # f = TF1("f", "x**x", 0, 1)
    # g = TF1("g", "x**x", 0, 1)



    # print(f.GetFormula())
    # print(g.GetFormula())
    
    # print(se.GetFormula())
    
    # ff = TF1("f", "SE_sc_sgn_centr", 0, 3)
    # ff.Draw()

    cf = TF1('f', getCorrRoot, 0, 3)

    cf.Draw()
    # cf = CorrelationHandler(num, den)

    # func = cf.eval

    # ff = TF1("ff", func, 0, 10)
    # ff.Draw()
    # # cf.corr_func.Draw()

    # # func(2)
    # input()

    




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('inFileName', metavar='text',
                        help='AnalysisResults.root file')
    parser.add_argument('oFileName', metavar='text',
                        help='AnalysisResults.root file')

    args = parser.parse_args()

    extract_raw_cf(inFileName=args.inFileName, oFileName=args.oFileName)
