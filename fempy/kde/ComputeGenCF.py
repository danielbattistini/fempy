'''
Script to create a tree dataset from binned dataset
'''
from curses import KEY_SAVE
import numpy as np
import argparse
from sys import exit

from utils.io import GetObjectFromFile
from utils.handle import GetFemtoDreamPairId

from ROOT import TFile, TTree, gRandom, TH1D, TKDE, RDataFrame, TF1, SetOwnership, gROOT, TCanvas

def GetCorrelationFunction(se, me, x):
    return se.Eval(x) / me.Eval(x)

def ExtractGenCF(inFileName, oFileName):
    """
    Extract the genuine correlation function from the raw CF, sidebands, minijets etc. CFs.

    Parameters
    ----------
    inFileName : str
        path of input datafile

    Returns
    -------
    int
        Description of return value

    """

    pairCombs = GetFemtoDreamPairId('sc')
    regions = {'sgn': '', 'sbl': 'SBLeft_', 'sbr': 'SBRight_'}
    kdeEstimates = {}

    oFile = TFile(oFileName, 'recreate')
    for pairKey, pairVal in pairCombs.items():
        for regKey, regVal in regions.items():
            rdfSE = RDataFrame(f'treeDpi_SE_{regKey}_{pairVal}', inFileName)
            dataSE = rdfSE.AsNumpy(['kStar'])
            
            kStarSE = dataSE['kStar']
            nPairsSE = len(kStarSE)

            kdeSE = TKDE(nPairsSE, kStarSE, 0.0, 3, "KernelType:Gaussian;Iteration:Adaptive;Mirror:MirrorAsymBoth;Binning:RelaxedBinning")

            kdeEstimates[f'SE_{pairKey}_{regKey}_centr'] = kdeSE.GetFunction().Clone()
            kdeEstimates[f'SE_{pairKey}_{regKey}_centr'].AddToGlobalList()
            kdeEstimates[f'SE_{pairKey}_{regKey}_centr'].SetName(f'SE_{pairKey}_{regKey}_centr')
            kdeEstimates[f'SE_{pairKey}_{regKey}_centr'].Write(f'SE_{pairKey}_{regKey}_centr')
            SetOwnership(kdeEstimates[f'SE_{pairKey}_{regKey}_centr'], False)

            kdeEstimates[f'SE_{pairKey}_{regKey}_upper'] = kdeSE.GetUpperFunction(0.683)
            kdeEstimates[f'SE_{pairKey}_{regKey}_upper'].AddToGlobalList()
            kdeEstimates[f'SE_{pairKey}_{regKey}_upper'].SetName(f'SE_{pairKey}_{regKey}_upper')
            kdeEstimates[f'SE_{pairKey}_{regKey}_upper'].Write(f'SE_{pairKey}_{regKey}_upper')
            SetOwnership(kdeEstimates[f'SE_{pairKey}_{regKey}_upper'], False)

            kdeEstimates[f'SE_{pairKey}_{regKey}_lower'] = kdeSE.GetLowerFunction(0.683)
            kdeEstimates[f'SE_{pairKey}_{regKey}_lower'].AddToGlobalList()
            kdeEstimates[f'SE_{pairKey}_{regKey}_lower'].SetName(f'SE_{pairKey}_{regKey}_lower')
            kdeEstimates[f'SE_{pairKey}_{regKey}_lower'].Write(f'SE_{pairKey}_{regKey}_lower')
            SetOwnership(kdeEstimates[f'SE_{pairKey}_{regKey}_lower'], False)
    
            # Mixed events
            rdfME = RDataFrame(f'treeDpi_ME_{regKey}_{pairVal}', inFileName)
            dataME = rdfME.AsNumpy(['kStar'])
            
            kStarME = dataME['kStar']
            nPairsME = len(kStarME)

            kdeME = TKDE(nPairsME, kStarME, 0.0, 3, "KernelType:Gaussian;Iteration:Adaptive;Mirror:MirrorAsymBoth;Binning:RelaxedBinning")

            kdeEstimates[f'ME_{pairKey}_{regKey}_centr'] = kdeME.GetFunction().Clone()
            kdeEstimates[f'ME_{pairKey}_{regKey}_centr'].AddToGlobalList()
            kdeEstimates[f'ME_{pairKey}_{regKey}_centr'].SetName(f'ME_{pairKey}_{regKey}_centr')
            kdeEstimates[f'ME_{pairKey}_{regKey}_centr'].Write(f'ME_{pairKey}_{regKey}_centr')
            SetOwnership(kdeEstimates[f'ME_{pairKey}_{regKey}_centr'], False)

            kdeEstimates[f'ME_{pairKey}_{regKey}_upper'] = kdeME.GetUpperFunction(0.683)
            kdeEstimates[f'ME_{pairKey}_{regKey}_upper'].AddToGlobalList()
            kdeEstimates[f'ME_{pairKey}_{regKey}_upper'].SetName(f'ME_{pairKey}_{regKey}_upper')
            kdeEstimates[f'ME_{pairKey}_{regKey}_upper'].Write(f'ME_{pairKey}_{regKey}_upper')
            SetOwnership(kdeEstimates[f'ME_{pairKey}_{regKey}_upper'], False)

            kdeEstimates[f'ME_{pairKey}_{regKey}_lower'] = kdeME.GetLowerFunction(0.683)
            kdeEstimates[f'ME_{pairKey}_{regKey}_lower'].AddToGlobalList()
            kdeEstimates[f'ME_{pairKey}_{regKey}_lower'].SetName(f'ME_{pairKey}_{regKey}_lower')
            kdeEstimates[f'ME_{pairKey}_{regKey}_lower'].Write(f'ME_{pairKey}_{regKey}_lower')
            SetOwnership(kdeEstimates[f'ME_{pairKey}_{regKey}_lower'], False)

    weightLeft = 0.51
    purity = 0.71

    def GetGenCorrelationFunction(distr, x, pairKey, regKey):
        se_sgn = distr[f'SE_{pairKey}_sgn_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sgn_{regKey}'].Integral(0, 3, 1e-3)
        me_sgn = distr[f'ME_{pairKey}_sgn_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sgn_{regKey}'].Integral(0, 3, 1e-3)

        se_sbl = distr[f'SE_{pairKey}_sbl_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sbl_{regKey}'].Integral(0, 3, 1e-3)
        me_sbl = distr[f'ME_{pairKey}_sbl_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sbl_{regKey}'].Integral(0, 3, 1e-3)

        se_sbr = distr[f'SE_{pairKey}_sbr_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sbr_{regKey}'].Integral(0, 3, 1e-3)
        me_sbr = distr[f'ME_{pairKey}_sbr_{regKey}'].Eval(x) / distr[f'SE_{pairKey}_sbr_{regKey}'].Integral(0, 3, 1e-3)

        cf_sgn = se_sgn/me_sgn
        cf_sbl = se_sbl/me_sbl
        cf_sbr = se_sbr/me_sbr

        cf_sb = weightLeft * cf_sbl + (1-weightLeft) * cf_sbr
        
        cf = (cf_sgn * purity - (1-purity)* cf_sb) /purity

        return cf_sgn


    cRaw = TCanvas("cRaw", "cRaw"regions
    fCFRawCentr = TF1("fCFRawCentr", cfRawCentr, 0, 3) 
    fCFRawUpper = TF1("fCFRawUpper", cfRawUpper, 0, 3) 
    fCFRawLower = TF1("fCFRawLower", cfRawLower, 0, 3) 

    fCFRawCentr.Draw()
    fCFRawUpper.Draw('same')
    fCFRawLower.Draw('same')
    cRaw.Write()
   
    # cSBL = TCanvas("cSBL", "cSBL", 600, 600)
    # cfSBLCentr = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbl_centr'], kdeEstimates[f'ME_pp_sbl_centr'], x[0])
    # cfSBLUpper = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbl_upper'], kdeEstimates[f'ME_pp_sbl_upper'], x[0])
    # cfSBLLower = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbl_lower'], kdeEstimates[f'ME_pp_sbl_lower'], x[0])
    # fCFSBLCentr = TF1("fCFSBLCentr", cfSBLCentr, 0, 3) 
    # fCFSBLUpper = TF1("fCFSBLUpper", cfSBLUpper, 0, 3) 
    # fCFSBLLower = TF1("fCFSBLLower", cfSBLLower, 0, 3) 

    # fCFSBLCentr.Draw()
    # fCFSBLUpper.Draw('same')
    # fCFSBLLower.Draw('same')
    # cSBL.Write()

    # cSBR = TCanvas("cSBR", "cSBR", 600, 600)
    # cfSBRCentr = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbr_centr'], kdeEstimates[f'ME_pp_sbr_centr'], x[0])
    # cfSBRUpper = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbr_upper'], kdeEstimates[f'ME_pp_sbr_upper'], x[0])
    # cfSBRLower = lambda x, par: GetGenCorrelationFunction(kdeEstimates[f'SE_pp_sbr_lower'], kdeEstimates[f'ME_pp_sbr_lower'], x[0])
    # fCFSBRCentr = TF1("fCFSBRCentr", cfSBRCentr, 0, 3) 
    # fCFSBRUpper = TF1("fCFSBRUpper", cfSBRUpper, 0, 3) 
    # fCFSBRLower = TF1("fCFSBRLower", cfSBRLower, 0, 3) 

    # fCFSBRCentr.Draw()
    # fCFSBRUpper.Draw('same')
    # fCFSBRLower.Draw('same')
    # cSBR.Write()

    
    
    
        # cfLambdaCentr = lambda x, par: GetGenCorrelationFunction(kdeEstimates, x[0], 'pp', 'centr')
    # cfLambdaUpper = lambda x, par: GetGenCorrelationFunction(kdeEstimates, x[0], 'pp', 'upper')
    # cfLambdaLower = lambda x, par: GetGenCorrelationFunction(kdeEstimates, x[0], 'pp', 'lower')

    # cfCentr = TF1("f", cfLambdaCentr, 0, 3) 
    # cfUpper = TF1("f", cfLambdaUpper, 0, 3) 
    # cfLower = TF1("f", cfLambdaLower, 0, 3) 
    
    # cfCentr.Draw()
    # cfUpper.Draw('same')
    # cfLower.Draw('same')

    cRaw.Write()

    oFile.Close()

    input()
    
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('inFileName', metavar='text',
                        help='AnalysisResults.root file')
    parser.add_argument('oFileName', metavar='text',
                        help='AnalysisResults.root file')

    args = parser.parse_args()

    ExtractGenCF(inFileName=args.inFileName, oFileName=args.oFileName)
