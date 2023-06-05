from ROOT import TFile, kBlue, kBlack, kRed, TCanvas, TLegend, kGreen

def GetNormFactor(se, me, fromVal, toVal):
    firstBin = se.FindBin(fromVal*0.9999)
    lastBin = se.FindBin(toVal*0.9999)
    
    return me.Integral(firstBin, lastBin) / se.Integral(firstBin, lastBin)


inFile = TFile('/home/daniel/an/DstarPi/20_luuksel/distr/Distr_data_nopc_kStarBW15MeV.root')
hME = {}
hSE = {}
hCF = {}
for region in ['sgn', 'sbr']:
    hSE[region] = inFile.Get(f'sc/SE/{region}/hCharmMassVsKStar0').ProjectionX(f'hSE_{region}')
    hME[region] = inFile.Get(f'sc/ME/{region}/hCharmMassVsKStar0').ProjectionX(f'hME_{region}')
    hSE[region].Sumw2()
    hSE[region].Rebin(4)
    hME[region].Rebin(4)
    norm = GetNormFactor(hSE[region], hME[region], 1500, 2000)
    print(norm)
    hCF[region] = hSE[region]/hME[region]
    hCF[region].Scale(norm)

# SB corr
cCompare = TCanvas('', '', 600, 600)
cCompare.DrawFrame(0, 0, 200, 5, ';k* (MeV/c);C(k*)')
purity = 0.78
hCF['sbcorr'] = 1/purity * ( hCF['sgn'] - (1-purity)*hCF['sbr'])
hCF['sbr'].SetLineColor(kBlue)
hCF['sbr'].SetLineWidth(2)
hCF['sbr'].Draw('same')

hCF['sgn'].SetLineColor(kRed)
hCF['sgn'].SetLineWidth(2)
hCF['sgn'].Draw('same')

hCF['sbcorr'].SetLineColor(kBlack)
hCF['sbcorr'].SetLineWidth(2)
hCF['sbcorr'].Draw("same")

sblessFile = TFile('/home/daniel/an/DstarPi/20_luuksel/cf/RawCF_data_nopc_kStarBW20MeV_gaus_charmMassBW0.2MeV_lowKStarCharmMassBW0.4MeV_until50MeV.root')
hGood = sblessFile.Get('sc/hCFPurityRewFromDataMinusBkg0')
hGood.SetLineColor(kGreen+2)
hGood.SetLineWidth(2)
hGood.Draw('same')

leg = TLegend(0.4, 0.75, 0.9, 0.9)
leg.AddEntry(hCF['sgn'], 'Signal')
leg.AddEntry(hCF['sbr'], 'SB right')
leg.AddEntry(hCF['sbcorr'], 'SB corrected')
leg.AddEntry(hGood, 'SB-less')
leg.Draw('same')

cCompare.SaveAs('/home/daniel/an/DstarPi/20_luuksel/comparisons/CompareSBApproach.pdf')