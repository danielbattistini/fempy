import argparse
import os
import yaml
from rich import print

from ROOT import TFile, TCanvas, TLegend, TLine, TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D, gStyle

import fempy
from fempy import logger as log
from fempy.utils.format import TranslateToLatex
from fempy.utils.io import Load
from fempy.utils import style

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfg')
args = parser.parse_args()

# Load configuration file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError:
        log.critical('Yaml file not loaded')

style.SetStyle()

for plot in cfg:
    plot = plot["plot"]

    panels = {'default': 1}
    if plot['ratio']['enable']:
        panels['ratio'] = len(panels)+1
    if plot['relunc']['enable']:
        panels['relunc'] = len(panels)+1

    # Load the objects to draw
    graphsStat = []
    graphsSyst = []
    legends = []
    drawOpts = []
    for inputCfg in plot["input"]:
        inFile = TFile(inputCfg['file'])

        inObjStat = Load(inFile, inputCfg['histostat'])
        inObjSyst = Load(inFile, inputCfg['histosyst'])

        if isinstance(inObjStat, TH1):
            inObjStat.SetDirectory(0)
            inObjStat.Rebin(inputCfg['rebin'])

            if inputCfg['normalize']:
                inObjStat.Scale(1./inObjStat.Integral())
            if inputCfg['normalizecf']:
                inObjStat.Scale(inputCfg['normalizecf'])
        
        statErr = TGraphAsymmErrors(inObjStat)
        systErr = TGraphAsymmErrors(inObjSyst)

        statErr.SetLineColor(style.GetColor(inputCfg['color']))
        statErr.SetMarkerColor(style.GetColor(inputCfg['color']))
        statErr.SetLineWidth(inputCfg.get('thickness', 1))
        statErr.SetMarkerStyle(inputCfg['markerstyle'])
        statErr.SetMarkerSize(inputCfg['markersize'])
        systErr.SetLineColor(style.GetColor(inputCfg['color']))
        systErr.SetMarkerColor(style.GetColor(inputCfg['color']))
        systErr.SetLineWidth(inputCfg.get('thickness', 1))
        systErr.SetMarkerStyle(inputCfg['markerstyle'])
        systErr.SetMarkerSize(inputCfg['markersize'])
        systErr.SetFillColor(style.GetColor(inputCfg['color']))
        if(inputCfg['fillstylesyst'] is not None):
            systErr.SetFillStyle(inputCfg['fillstylesyst'])
            systErr.SetFillColorAlpha(style.GetColor(inputCfg['color']), inputCfg['fillalphasyst'])
        
        for iPoint in range(inObjStat.GetNbinsX()):
            errX = inObjStat.GetBinWidth(iPoint)/4
            statErr.SetPointEXlow(iPoint, errX*2)
            statErr.SetPointEXhigh(iPoint, errX*2)
            systErr.SetPointEXlow(iPoint, errX)
            systErr.SetPointEXhigh(iPoint, errX)        
        
        graphsStat.append(statErr)
        graphsSyst.append(systErr)
        
        legends.append(inputCfg['legend'])

        drawOpts.append(inputCfg.get('drawopt', 'p' if isinstance(statErr, TH1) else 'pe'))
    # Define the canvas
    nPanelsX, nPanelsY = fempy.utils.GetNPanels(len(panels))
    cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
    cPlot.Divide(nPanelsX, nPanelsY)
    pad = cPlot.cd(1)
    pad.SetLogx(plot["opt"]["logx"])
    pad.SetLogy(plot["opt"]["logy"])
    if(plot["opt"]["padtopmargin"]):
        pad.SetTopMargin(plot["opt"]["padtopmargin"])
    if(plot["opt"]["padbottommargin"]):
        pad.SetBottomMargin(plot["opt"]["padbottommargin"])
    if(plot["opt"]["padrightmargin"]):
        pad.SetRightMargin(plot["opt"]["padrightmargin"])
    if(plot["opt"]["padleftmargin"]):
        pad.SetLeftMargin(plot["opt"]["padleftmargin"])
    if(plot['opt']['ytitleoffset'] is not None):
        print(plot['opt']['ytitleoffset'])
        gStyle.SetTitleOffset(plot['opt']['ytitleoffset'],"Y")

    fx1 = plot['opt']['rangex'][0]
    fy1 = plot['opt']['rangey'][0]
    fx2 = plot['opt']['rangex'][1]
    fy2 = plot['opt']['rangey'][1]
    pad.DrawFrame(fx1, fy1, fx2, fy2, TranslateToLatex(plot['opt']['title']))

    legx1 = plot['opt']['leg']['posx'][0]
    legy1 = plot['opt']['leg']['posy'][0]
    legx2 = plot['opt']['leg']['posx'][1]
    legy2 = plot['opt']['leg']['posy'][1]
    leg = TLegend(legx1, legy1, legx2, legy2)

    for iObj, (inObjStat, inObjSyst, legend) in enumerate(zip(graphsStat, graphsSyst, legends)):
        inObjStat.Draw(inputCfg['drawoptstat'])
        inObjSyst.Draw(inputCfg['drawoptsyst'])

        # Compute statistics for hist in the displayed range
        if isinstance(inObjStat, TH1):
            firstBin = inObjStat.FindBin(plot['opt']['rangex'][0]*1.0001)
            lastBin = inObjStat.FindBin(plot['opt']['rangex'][1]*0.9999)
            inObjStat.GetXaxis().SetRange(firstBin, lastBin)
            print(f'{legend}: mean = {inObjStat.GetMean()} sigma = {inObjStat.GetStdDev()}')
            if plot['opt']['leg']['mean']:
                legend += f';  #mu={inObjStat.GetMean():.3f}'
            if plot['opt']['leg']['sigma']:
                legend += f';  #sigma={inObjStat.GetStdDev():.3f}'
        leg.AddEntry(inObjStat, legend, 'lp')
        
    if(plot['opt']['leg']['center']): 
        leg.SetHeader(TranslateToLatex(plot['opt']['leg']['header']), 'C')
    else: 
        leg.SetHeader(TranslateToLatex(plot['opt']['leg']['header']))

    leg.Draw()

    cPlot.Modified()
    cPlot.Update()

    # save canvas
    for ext in plot["opt"]["ext"]:
        cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')































#import argparse
#import os
#import yaml
#import numpy as np
#from rich import print
#
#from ROOT import TFile, TCanvas, TH1, TLegend, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D, TBox, TLatex, TText
#
#import fempy
#from fempy import logger as log
#from fempy.utils.format import TranslateToLatex
#from fempy.utils.io import Load
#from fempy.utils import style
#
#parser = argparse.ArgumentParser(description='Arguments')
#parser.add_argument('cfg')
#args = parser.parse_args()
#
## Load configuration file
#with open(args.cfg, "r") as stream:
#    try:
#        cfg = yaml.safe_load(stream)
#    except yaml.YAMLError:
#        log.critical('Yaml file not loaded')
#
#style.SetStyle()
#
#for plot in cfg:
#    plot = plot["plot"]
#
#    nInputs = len(plot['input']) 
#    panels = {'default': 1}
#    if(nInputs>1):
#        panels['uncratio'] = len(panels)+1
#
#    # Define the canvas
#    nPanelsX, nPanelsY = fempy.utils.GetNPanels(len(panels))
#    cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
#    cPlot.Divide(nPanelsX, nPanelsY)
#    pad = cPlot.cd(1)
#    pad.SetLogx(plot["opt"]["logx"])
#    pad.SetLogy(plot["opt"]["logy"])
#
#    fx1 = plot['opt']['rangex'][0]
#    fy1 = plot['opt']['rangey'][0]
#    fx2 = plot['opt']['rangex'][1]
#    fy2 = plot['opt']['rangey'][1]
#    print(fx1)
#    pad.DrawFrame(fx1, fy1, fx2, fy2, TranslateToLatex(plot['opt']['title']))   
#
#    boxes = []
#    inObjStatGraphs = []
#    legends = []
#    
#    # Load the objects to draw
#    for inputCfg in plot["input"]:
#        inFile = TFile(inputCfg['file'])
#
#        inObjStat = Load(inFile, inputCfg['name'])
#        legends.append(inputCfg['legend'])
#        print(inputCfg['legend'])
#        if(isinstance(inObjStat, TH1)):
#            inObjStat.SetDirectory(0)
#            inObjStat.Rebin(inputCfg['rebin'])        
#            if inputCfg['normalize']:
#                inObjStat.Scale(1./inObjStat.Integral())
#
#        binCenters = []
#        binWidths = []
#        binContents = []
#        binErrors = []
#        for iBin in range(inObjStat.GetNbinsX()):
#            binCenters.append(inObjStat.GetBinCenter(iBin + 1))
#            binWidths.append(inObjStat.GetBinWidth(iBin + 1))
#            binContents.append(inObjStat.GetBinContent(iBin +1))
#            binErrors.append(inObjStat.GetBinError(iBin + 1))
#            if(nInputs==1): 
#                if(inputCfg['errbarfillcolor'] is not None):
#                    box = TBox(binCenters[-1] - binWidths[-1]/4, binContents[-1] - binErrors[-1], 
#                               binCenters[-1] + binWidths[-1]/4, binContents[-1] + binErrors[-1])
#                    box.SetFillStyle(inputCfg['errbarfillstyle'])
#                    box.SetFillColorAlpha(inputCfg['errbarfillcolor'], inputCfg['errbarfillalpha'])
#                    #box.Draw('same')
#                    boxes.append(box)
#        binCentersNp = np.array(binCenters, dtype='double')
#        binContentsNp = np.array(binContents, dtype='double')
#        binXErrLowNp = np.array([binWidth/2 for binWidth in binWidths], dtype='double')
#        binXErrUppNp = np.array([binWidth/2 for binWidth in binWidths], dtype='double')
#        binYErrLowNp = np.array([binErr/2 for binErr in binErrors], dtype='double')
#        binYErrUppNp = np.array([binErr/2 for binErr in binErrors], dtype='double')
#        inObjStatGraph = TGraphAsymmErrors(inObjStat.GetNbinsX(), binCentersNp, binContentsNp, binXErrLowNp, 
#                                       binXErrUppNp, binYErrLowNp, binYErrUppNp)   
#        inObjStatGraph.GetXaxis().SetLimits(fx1, fx2) 
#        inObjStatGraph.GetYaxis().SetLimits(fy1, fy2) 
#        inObjStatGraph.GetXaxis().SetRangeUser(fx1, fx2) 
#        inObjStatGraph.GetYaxis().SetRangeUser(fy1, fy2) 
#        inObjStatGraph.SetLineColor(style.GetColor(inputCfg['color']))
#        inObjStatGraph.SetMarkerColor(style.GetColor(inputCfg['color']))
#        inObjStatGraph.SetMarkerStyle(inputCfg['markerstyle'])
#        inObjStatGraph.SetMarkerSize(inputCfg['markersize'])
#        inObjStatGraph.SetTitle(TranslateToLatex(plot['opt']['title']))
#        inObjStatGraphs.append(inObjStatGraph)  
#        print('ciao')
#        inObjStatGraph.Draw('same ap')
#    
#        cPlot.Modified()
#        cPlot.Update()
#    
#    if(nInputs==1): 
#        for iBox in boxes:
#            cPlot.Modified()
#            cPlot.Update()
#            iBox.Draw('same')
#        
#    #print(len(inObjStatGraphs))
#    #for inObjStatGraph in inObjStatGraphs:
#    #    cPlot.Modified()
#    #    cPlot.Update()
#    #    inObjStatGraph.Draw('same ap')
#
#    cPlot.Modified()
#    cPlot.Update()
#        
#    legx1 = plot['opt']['leg']['posx'][0]
#    legy1 = plot['opt']['leg']['posy'][0]
#    legx2 = plot['opt']['leg']['posx'][1]
#    legy2 = plot['opt']['leg']['posy'][1]
#    leg = TLegend(legx1, legy1, legx2, legy2)
#        
#    print(inObjStatGraphs)
#    print(legends)
#    print(len(boxes))
#    for (inObjStatGraph, legend) in zip(inObjStatGraphs, legends):
#        leg.AddEntry(inObjStatGraph, legend, 'lp')    
#    
#    leg.SetHeader(TranslateToLatex(plot['opt']['leg']['header']), 'C')
#    if(plot['opt']['leg']['enable']):
#        leg.Draw()
#    
#    for text in plot['opt']['description']:
#        if('#' in text['text']):
#            tl = TLatex()
#            tl.SetTextSize(text['textsize'])
#            tl.SetTextFont(text['textfont'])
#            tl.DrawLatexNDC(text['position'][0], text['position'][1], text['text'])
#            cPlot.Modified()
#            cPlot.Update()
#        else:
#            tl = TLatex()
#            tl.SetTextSize(text['textsize'])
#            tl.SetTextFont(text['textfont'])
#            tl.DrawTextNDC(text['position'][0], text['position'][1], text['text'])
#            cPlot.Modified()
#            cPlot.Update()
#        
#    if(nInputs==1):
#        pad = cPlot.cd(panels['uncratio'])
#        pad.SetLogx(plot['ratio']['logx'])
#        pad.SetLogy(plot['ratio']['logy'])
#        x1 = plot['opt']['rangex'][0]
#        y1 = plot['ratio']['rangey'][0]
#        x2 = plot['opt']['rangex'][1]
#        y2 = plot['ratio']['rangey'][1]
#        frame = pad.DrawFrame(x1, y1, x2, y2, 'Rel. Unc. comparison')
#        frame.GetYaxis().SetTitle('Ratio')
#        hDen = inObjStats[0].Clone()
#        hDen.Rebin(plot['ratio']['rebin'])
#        hDen.Sumw2()
#    
#        if isinstance(inObjStat, TH1):
#            for inObjStat in inObjStats[1:]:
#                hRatio = inObjStat.Clone()
#                hRatio.Rebin(plot['ratio']['rebin'])
#                hRatio.Divide(hDen)
#                hRatio.Draw('same pe')
#        else:
#            log.error('Ratio for type %s is not implemented. Skipping this object', type(inObjStat))
#            continue
#    
#        line = TLine(plot['opt']['rangex'][0], 1, plot['opt']['rangex'][1], 1)
#        line.SetLineColor(13)
#        line.SetLineStyle(7)
#        line.Draw('same pe')
#    
#    cPlot.Modified()
#    cPlot.Update()
#
#    # save canvas
#    for ext in plot["opt"]["ext"]:
#        cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')