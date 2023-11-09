from ROOT import gRandom, TH1D, TCanvas

gRandom.SetSeed(42)

nIter = 2000000
hRatio=TH1D('hRatio', 'hRatio', 200, 0, 10)


c1 = TCanvas("c1")
muNum = 11
muDen = 2
for iIter in range(nIter):
    num = gRandom.Poisson(muNum)
    den = gRandom.Poisson(muDen)
    if den > 0:
        hRatio.Fill(num/den)
hRatio.Draw()

c2 = TCanvas("c2")
hGauss = TH1D('hGauss', 'hGauss', 22, 0, 22)
hPoisson = TH1D('hPoisson', 'hPoisson', 22, 0, 22)

mu = 11
for iIter in range(nIter):
    hPoisson.Fill(gRandom.Poisson(mu))
    hGauss.Fill(gRandom.Gaus(mu, mu**0.5))
hGauss.Draw()
hPoisson.Draw('same')
