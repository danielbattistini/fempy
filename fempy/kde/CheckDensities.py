import sys

from ROOT import TFile, TCanvas, TF1, kBlack, kBlue, kRed

inFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/densities/Densities_test_kern-Gaussian_iter-Adaptive_mirror-NoMirror_binning-RelaxedBinning_noraremultbins.root"
inFile = TFile(inFileName)
kStarMax = 3000

# get n ofg mult bins
nMultBins = 25


# get densities for each mult bin
pId = 'sc'
rId = 'sgn'
c = TCanvas("cSlices", "cSlices", 600, 600)
inFile.ls()

fList = []
# for eMultBin, iMultBin in enumerate([5, 10, 15, 20]):
for eMultBin, iMultBin in enumerate(range(25)):
    multMin = iMultBin
    multMax = iMultBin+1
    
    f = inFile.Get(f'fME_mult_{multMin:.0f}_{multMax:.0f}_{pId}_{rId}_centr')
    if f is None:
        continue
    f.SetLineColor(kBlack)
    if eMultBin == 0:
        f.GetYaxis().SetRangeUser(0, 0.0009)
        f.Draw()
    else:
        f.Draw('same')

    fList.append(f)

hWeights = inFile.Get(f'hWeights_{pId}_{rId}')
norm = sum([hWeights.GetBinContent(iSlice) * fList[iSlice].Integral(0, kStarMax, 1e-4) for iSlice in range(nMultBins)])

def reweight(x, par):
    return sum([hWeights.GetBinContent(iSlice) * fList[iSlice].Eval(x[0]) for iSlice in range(nMultBins)])/norm

cfMERew = TF1("rew", reweight, 0, kStarMax)
cWeights = TCanvas("cWeights", "", 600, 600)
cfMERew.Draw()

fME = inFile.Get(f'fME_{pId}_{rId}_centr')
fME.SetLineColor(kBlue)
fME.Draw("same")

hWeights = inFile.Get(f'hWeights_{pId}_{rId}')
hMEProjNorm = inFile.Get(f'hMEProjNorm')
hSEProjNorm = inFile.Get(f'hSEProjNorm')

cMult = TCanvas("cMult", "", 600, 600)
hWeights.SetLineColor(kBlack)
hWeights.Draw()

hMEProjNorm.SetLineColor(kBlue)
hMEProjNorm.Scale(4)
hMEProjNorm.Draw("same hist")

hSEProjNorm.SetLineColor(kRed)
hSEProjNorm.Scale(4)
hSEProjNorm.Draw("same hist")


