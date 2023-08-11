'''test
'''
from ROOT import TF1, TFile

from utils.correlation_handler import CorrelationHandler
from utils.io import GetObjectFromFile

se = TF1('se', 'x+x*x', 0, 3)
me = TF1('me', 'x+0.9*x*x', 0, 3)

cf = CorrelationHandler(se, me)
mycf = cf.corr_func
mycf.Draw()

inFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/GenCF_red.root'
inFile = TFile(inFileName)

se = GetObjectFromFile(inFile, 'SE_pp_sgn_centr')
me = GetObjectFromFile(inFile, 'ME_pp_sgn_centr')
