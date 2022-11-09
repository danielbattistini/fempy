import argparse
import sys
import os
import yaml

from ROOT import TFile, TCanvas, TLegend, TLine

from utils.io import GetObjectFromFile
from utils.format import TranslateToLatex, FigInit

FigInit()

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfg')
parser.add_argument('oFileName')
args = parser.parse_args()

with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

doRatio = cfg['ratio']['enable']

if doRatio:
    cPlot = TCanvas("cPlot", "cPlot", 1200, 600)
    cPlot.Divide(2, 1)
    cPlot.cd(1)
else:
    cPlot = TCanvas("cPlot", "cPlot", 600, 600)

inObjs = []
legends = []
for input in cfg["input"]:
    inFile = TFile(input['file'])

    inObj = GetObjectFromFile(inFile, input['name'])
    if inObj == None:
        print(f'Error: cannot load {inFile}:{input["name"]}. Exit!')
        sys.exit()

    inObj.Rebin(cfg['opt']['rebin'])
    inObj.SetDirectory(0)
    inObjs.append(inObj)
    legends.append(input['legend'])

legx1 = cfg['opt']['leg']['posx'][0]
legy1 = cfg['opt']['leg']['posy'][0]
legx2 = cfg['opt']['leg']['posx'][1]
legy2 = cfg['opt']['leg']['posy'][1]
leg = TLegend(legx1, legy1, legx2, legy2)

for iObj, (inObj, legend) in enumerate(zip(inObjs, legends)):
    inObj.SetStats(0)
    inObj.SetLineWidth(2)
    if iObj == 0:
        inObj.SetTitle(cfg['opt']['title'])
        inObj.GetXaxis().SetRangeUser(
            cfg['opt']['rangex'][0], cfg['opt']['rangex'][1])
        inObj.GetYaxis().SetRangeUser(
            cfg['opt']['rangey'][0], cfg['opt']['rangey'][1])
        inObj.Draw('plc pmc')
    else:
        inObj.Draw("same plc pmc")
    leg.AddEntry(inObj, legend)
leg.SetTextSize(0.03)
leg.SetTextSize(0.03)

leg.SetHeader(TranslateToLatex(cfg['opt']['leg']['header']), 'C')
leg.Draw()



# ratio
if doRatio:
    cPlot.cd(2)

    hDen = inObjs[0].Clone()
    hDen.Rebin(cfg['ratio']['rebin'])
    hRatio = inObjs[1].Clone()
    hRatio.Rebin(cfg['ratio']['rebin'])
    
    hRatio.Sumw2()
    
    hRatio.Divide(hDen)


    hRatio.SetTitle(cfg['opt']['title'])
    hRatio.GetYaxis().SetTitle('Ratio')
    hRatio.GetYaxis().SetRangeUser(cfg['ratio']['rangey'][0], cfg['ratio']['rangey'][1])
    hRatio.GetXaxis().SetRangeUser(cfg['opt']['rangex'][0], cfg['opt']['rangex'][1])
    hRatio.Draw()

    line = TLine(cfg['opt']['rangex'][0], 1, cfg['opt']['rangex'][1], 1)
    line.SetLineColor(13)
    line.SetLineStyle(9)

    line.Draw('same')

for ext in cfg["opt"]["ext"]:
    cPlot.SaveAs(f'{os.path.splitext(args.oFileName)[0]}.{ext}')
