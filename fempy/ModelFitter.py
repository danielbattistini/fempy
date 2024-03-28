'''
Script to perform the fit on a correlation function.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 ModelFitter.py cfg.yml

'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, TH1D
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/ModelFitter.hxx"')
from ROOT import ModelFitter

from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import ChangeUnits

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--debug', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Load input file with data and mc CF
inFileFit = TFile(cfg['infile'])

# Define the output file
oFileBaseName = cfg['ofilebasename']
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

fileLines = []
fitCfgIdxs = []
with open(args.cfg, 'r') as file:
    for lineNumber, line in enumerate(file, start=1):
        if('cfpath' in line):
            fitCfgIdxs.append(lineNumber)
        firstNonBlankChar = next((char for char in line if not char.isspace()), None)
        if(firstNonBlankChar == '#'): continue
        fileLines.append(line.rstrip())  # Strip to remove leading/trailing whitespaces

oFileNameCfg = os.path.join(cfg['odir'], oFileBaseName + '_cfg.txt')
with open(oFileNameCfg, 'w') as file:
    for line in fileLines:
        file.write(line)
        file.write('\n')

modelFitters = []

# for loop over the correlation functions
for fitcf in cfg['fitcfs']:
    
    # change unity of measure of histograms from GeV to MeV
    fitHisto = ChangeUnits(Load(inFileFit, fitcf['cfpath']), 1000)

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    if('rejectrange' in fitcf):
        lowRejectRange = fitcf['rejectrange'][0]
        uppRejectRange = fitcf['rejectrange'][1]
        modelFitters.append(ModelFitter(fitHisto, lowFitRange, uppFitRange, lowRejectRange, uppRejectRange))
    else: 
        modelFitters.append(ModelFitter(fitHisto, lowFitRange, uppFitRange))
        
    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    compsToFile = []
    # for loop over the functions entering in the model
    for iFunc, func in enumerate(fitcf['model']):
        
        if('isbaseline' in func):
            if(func['isbaseline']):
                modelFitters[-1].SetBaselineIdx(iFunc)
        if('savetofile' in func):
            compsToFile.append(iFunc)
        
        # fit function parameters initialization
        initPars = []
    
        if('splinehisto' in func['funcname']):
            histoFile = TFile(func['histofile'])
            splinedHisto = ChangeUnits(Load(histoFile, func['histopath']), 1000)
            if('rebin' in func):
                splinedHisto.Rebin(func['rebin'])
            if('p0' in func):
                initPars = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3]),
                            (func['p0'][0], func['p0'][1], func['p0'][2], func['p0'][3])]
            else:
                initPars = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
            modelFitters[-1].AddSplineHisto(func['funcname'], splinedHisto, initPars, func['addmode'], func['onbaseline'])
            cSplinedHisto = TCanvas(f'c{func["funcname"]}', '', 600, 600)
            modelFitters[-1].DrawSpline(cSplinedHisto, splinedHisto)
            oFile.cd(fitcf['fitname'])
            cSplinedHisto.Write()
            
        #elif('TF1' in func['funcname']):
        #    funcFile = TFile(func['TF1file'])
        #    fixedTF1 = Load(funcFile, func['TF1path'])
        #    initPars = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
        #    modelFitters[-1].AddTF1(func['funcname'], fixedTF1, initPars, func['addmode'], func['onbaseline'])
        #    cFixedTF1 = TCanvas(f'c{func["funcname"]}', '', 600, 600)
        #    modelFitters[-1].DrawTF1Comp(cFixedTF1, fixedTF1)
        #    oFile.cd(fitcf['fitname'])
        #    cFixedTF1.Write()
        
        else:  
            #if('spline3' in func['funcname']):
            #    for nKnot, xKnot in enumerate(func['xknots']):
            #        initPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
            #    for nKnot, xKnot in enumerate(func['xknots']):
            #        nBin = prefitHisto.FindBin(xKnot)
            #        yKnot = prefitHisto.GetBinContent(nBin)
            #        initPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*30, yKnot + (yKnot/100)*30])
            #else:
            if('fixparsfromfuncts' in func):
                histoFuncFile = TFile(func['histofuncfile'])
                histoFuncParams = Load(histoFuncFile, func['histofuncpath'])
                for iPar in func['fixparsfromfuncts']:
                    func[f'p{iPar[0]}'] = [f'p{iPar[0]}_prefit', 
                                           histoFuncParams.GetBinContent(iPar[1]), 0, -1]

            initPars = [(func[f'p{iPar}'][0], func[f'p{iPar}'][1], func[f'p{iPar}'][2], 
                         func[f'p{iPar}'][3]) for iPar in range(func['npars'])]

            if('lambdapar' in func):
                lambdaParam = [("lambdapar_" + func['funcname'], func['lambdapar'], 0, -1)]
                initPars = lambdaParam + initPars
                modelFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
            if('lambdagen' in func):
                lambdaGen = [("lambda_gen_" + func['funcname'], func['lambdagen'], 0, -1)]
                initPars = lambdaGen + initPars
                modelFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
            if('norm' in func):
                normParam = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
                initPars = normParam + initPars
                modelFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
                
    # perform the fit and save the result
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].Fit()
    cFit = TCanvas('cFit', '', 600, 600)
    if('drawsumcomps' in fitcf):
        modelFitters[-1].Draw(cFit, fitcf['drawsumcomps'])
    else:
        modelFitters[-1].Draw(cFit)
        
    modelFitters[-1].DrawLegend(cFit, fitcf['legcoords'][0], fitcf['legcoords'][1], fitcf['legcoords'][2], 
                             fitcf['legcoords'][3], fitcf['legentries'])
    cFit.Write()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    fitFunction.Write()
    for compToFile in compsToFile:
        if('spline' not in fitcf['model'][iFunc]['funcname']):
            modelFitters[-1].GetComponent(compToFile).Write(func['funcname'])
        modelFitters[-1].GetComponentPars(compToFile).Write('h' + fitcf['model'][compToFile]['funcname'][0].upper() + 
                                                            fitcf['model'][compToFile]['funcname'][1:])
    
    with open(oFileNameCfg, 'a') as file:
        file.write('-----------------------------------')
        file.write('\n')
        file.write('Parameters obtained from the fit')    
        file.write('\n')
        for iPar in range(fitFunction.GetNpar()):
            file.write(fitFunction.GetParName(iPar) + ": " + str(fitFunction.GetParameter(iPar)))
            file.write('\n')
    
    #pdfFileName = fitcf['fitname'] + cfg["suffix"] + ".pdf"
    #pdfFilePath = os.path.join(cfg['odir'], pdfFileName) 
    #cFit.SaveAs(pdfFilePath)

oFile.Close()
print(f'Config saved in {oFileNameCfg}')
print(f'output saved in {oFileName}')