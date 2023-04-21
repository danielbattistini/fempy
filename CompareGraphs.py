import argparse
import sys
import os
import yaml
from rich import print
# todo: implement colors

from ROOT import TFile, TCanvas, TLegend, TLine, TH1, TGraph, TGraphAsymmErrors, gROOT, kRed, kBlack, kAzure

from fempy.utils.io import GetObjectFromFile
from fempy.utils.format import TranslateToLatex, FigInit
import fempy

colors = {
    'kRed': kRed,
    'kBlack': kBlack,
    'kAzure': kAzure,
}

FigInit()

def GetCanvasSplitting(nPanels):
    if nPanels < 3:
        return (nPanels, 1)
    elif nPanels < 8:
        return (round(nPanels/2), 2)
    else:
        fempy.error("not implemented")


def CompareGraphs(cfgName):
    with open(cfgName, "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()


        # gROOT.SetBatch(False)
        for plot in cfg:
            plot = plot["plot"]

            panels = {'default': 1}
            if plot['ratio']['enable']:
                panels['ratio'] = len(panels)+1
            if plot['relunc']['enable']:
                panels['relunc'] = len(panels)+1

            nPanelsX, nPanelsY = fempy.utils.GetNPanels(len(panels))
            cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
            cPlot.Divide(nPanelsX, nPanelsY)
            pad = cPlot.cd(1)

            inObjs = []
            legends = []
            for inputCfg in plot["input"]:
                inFile = TFile(inputCfg['file'])

                inObj = GetObjectFromFile(inFile, inputCfg['name'])
                if inObj == None:
                    fempy.error(f'cannot load {inFile}:{inputCfg["name"]}')

                if isinstance(inObj, TH1):
                    inObj.SetDirectory(0)
                    inObj.Rebin(inputCfg['rebin'])
                    inObj.SetLineColor(colors[inputCfg['color']])
                    inObj.SetMarkerColor(colors[inputCfg['color']])

                    print(inObj.GetNbinsX())
                    if inputCfg['normalize']:
                        inObj.Scale(1./inObj.Integral())
                inObjs.append(inObj)
                legends.append(inputCfg['legend'])

            legx1 = plot['opt']['leg']['posx'][0]
            legy1 = plot['opt']['leg']['posy'][0]
            legx2 = plot['opt']['leg']['posx'][1]
            legy2 = plot['opt']['leg']['posy'][1]
            leg = TLegend(legx1, legy1, legx2, legy2)

            for iObj, (inObj, legend) in enumerate(zip(inObjs, legends)):
                inObj.SetStats(0)
                inObj.SetLineWidth(2)
                if iObj == 0:
                    inObj.SetTitle(TranslateToLatex(plot['opt']['title']))
                    inObj.GetXaxis().SetRangeUser(
                        plot['opt']['rangex'][0], plot['opt']['rangex'][1])
                    inObj.GetYaxis().SetRangeUser(
                        plot['opt']['rangey'][0], plot['opt']['rangey'][1])
                    if isinstance(inObj, TGraph):
                        inObj.Draw('ap')
                    else:
                        inObj.Draw('')
                else:
                    inObj.Draw("same")
                leg.AddEntry(inObj, legend)
            leg.SetTextSize(0.03)
            leg.SetTextSize(0.03)

            leg.SetHeader(TranslateToLatex(plot['opt']['leg']['header']), 'C')
            leg.Draw()

            pad.SetLogx(plot["opt"]["logx"])
            pad.SetLogy(plot["opt"]["logy"])

            # ratio
            if plot['ratio']['enable']:
                pad = cPlot.cd(panels['ratio'])

                hDen = inObjs[0].Clone()
                hDen.Rebin(plot['ratio']['rebin'])

                hRatios = []
                for iInObj, inObj in enumerate(inObjs[1:]):
                    hRatio = inObj.Clone()
                    hRatio.Rebin(plot['ratio']['rebin'])
                    hRatio.Sumw2()

                    hRatio.Divide(hDen)

                    hRatio.SetTitle(TranslateToLatex(plot['opt']['title']))
                    hRatio.GetYaxis().SetTitle('Ratio')
                    hRatio.GetYaxis().SetRangeUser(plot['ratio']['rangey'][0], plot['ratio']['rangey'][1])
                    hRatio.GetXaxis().SetRangeUser(plot['opt']['rangex'][0], plot['opt']['rangex'][1])

                    hRatios.append(hRatio)

                    if iInObj == 0:
                        hRatio.Draw('')
                    else:
                        hRatio.Draw(' same')

                line = TLine(plot['opt']['rangex'][0], 1, plot['opt']['rangex'][1], 1)
                line.SetLineColor(13)
                line.SetLineStyle(9)

                line.Draw('same')
            
            # relative uncertainties
            if plot['relunc']['enable']:
                pad = cPlot.cd(panels['relunc'])
                pad.SetLeftMargin(0.16)

                inObjsUncs = []
                for iObj, inObj in enumerate(inObjs):
                    if isinstance(inObj, TH1):
                        inObjUnc = inObj.Clone()
                        for iBin in range(inObj.GetNbinsX()):
                            print(inObj.GetBinError(iBin+1)/inObj.GetBinContent(iBin+1))
                            inObjUnc.SetBinContent(iBin+1, 100 * inObj.GetBinError(iBin+1)/inObj.GetBinContent(iBin+1))
                            inObjUnc.SetBinError(iBin+1, 0)

                        inObjsUncs.append(inObjUnc)
                    elif isinstance(inObj, TGraph) and False:
                        for iPoint in range(inObj.GetN()):

                            if isinstance(inObjUnc, TGraphAsymmErrors):
                                x = inObj.GetPointX(iPoint)
                                y = inObj.GetPointY(iPoint)
                                xUncUpper = inObj.GetErrorXhigh(iPoint)
                                xUncLower = inObj.GetErrorXlow(iPoint)
                                yUncUpper = inObj.GetErrorYhigh(iPoint)
                                yUncLower = inObj.GetErrorYlow(iPoint)

                                yAvgUnc = (yUncLower + yUncUpper) / 2
                                inObjUnc.SetPoint(iPoint, x,  100 * yAvgUnc/y)

                                inObjUnc.SetPointError(iPoint, xUncLower, xUncUpper, 0, 0)
                            elif isinstance(inObjUnc, TGraphAsymmErrors):
                                x = inObj.GetPointX(iPoint)
                                y = inObj.GetPointY(iPoint)
                                xUnc = inObj.GetErrorX(iPoint)
                                yUnc = inObj.GetErrorY(iPoint)

                                inObjUnc.SetPoint(iPoint, x,  1000 * yUnc/y)
                                inObjUnc.SetPointError(iPoint, xUnc, 0)
                            else:
                                fempy.error('gne')

                    if iObj == 0:
                        inObjUnc.GetXaxis().SetRangeUser(plot['opt']['rangex'][0], plot['opt']['rangex'][1])
                        inObjUnc.GetYaxis().SetRangeUser(plot['relunc']['rangey'][0], plot['relunc']['rangey'][1])
                        inObjUnc.GetYaxis().SetTitle('Relative uncertainty (%)')
                        if isinstance(inObjUnc, TGraph):
                            inObjUnc.Draw('pa')
                        else:
                            print("drawwww")
                            inObjUnc.Draw('')
                    else:
                        inObjUnc.Draw('same')
                    inObjsUncs[0].SetLineColor(2)
                    # inObjsUncs[0].Draw()

                    # cPlot.Update()
                    # cPlot.Modified()

            # save canvas
            for ext in plot["opt"]["ext"]:
                cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')
    input()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('cfg')
    args = parser.parse_args()


        
    CompareGraphs(args.cfg)
