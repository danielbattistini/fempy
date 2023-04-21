from ROOT import TFile, TCanvas, TLegend, kRed, kBlack, gROOT,  kBlue, kGreen, kAzure

def MJFactorizatonPlots(colors):
    gROOT.SetBatch(True)
    baseDir = '/home/daniel/an/mjfactorization/sim'
    cMJ = TCanvas('cMJ', '', 600, 600)
    cMJ.DrawFrame(0, 0, 3, 5, ';k* (GeV/c);C(k*)')
    leg = TLegend(0.2, 0.7, 0.9, 0.9)
    hCF = []
    for mode, color in colors.items():
        inFile = TFile(f'{baseDir}/AnalysisResults_{mode}.root')
        hSE = inFile.Get('hSE')
        hME = inFile.Get('hME')
        hSE.SetDirectory(0)
        hME.SetDirectory(0)
        inFile.Close()
        print(hSE)
        if hSE == None:
            hSE = inFile.Get('hPairSE_22_211')
        if hME == None:
            hME = inFile.Get('hPairME_22_211')
        
        hSE.Rebin(100)
        hME.Rebin(100)

        hCF.append(hSE.Clone(f'hCF_{mode}'))
        hCF[-1].Sumw2()
        hCF[-1].Divide(hME)
        norm = [1, 1.5]
        hCF[-1].Scale(hME.Integral(hME.FindBin(norm[0]*1.0001), hME.FindBin(norm[1]*0.9999))/hSE.Integral(hSE.FindBin(norm[0]*1.0001), hSE.FindBin(norm[1]*0.9999)))
        hCF[-1].SetLineColor(color)
        hCF[-1].SetLineWidth(2)
        leg.AddEntry(hCF[-1], mode)
        hCF[-1].Draw('same')

    print(hCF)
    leg.Draw('same')
    cMJ.Modified()
    cMJ.Update()
    cMJ.SaveAs(f'{baseDir}/AnalysisResults_{"_".join(colors.keys())}.png')
    cMJ.SaveAs(f'{baseDir}/AnalysisResults_{"_".join(colors.keys())}.pdf')