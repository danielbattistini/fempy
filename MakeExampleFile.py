from ROOT import TFile, TH1D, gRandom

oFile = TFile("AnalysisResultsExample.root", 'recreate')
oFile.mkdir('very/complicated/path/')
oFile.cd('very/complicated/path/')
h1 = TH1D('h1', 'h1', 100, 0, 10)
h2 = TH1D('h2', 'h2', 100, 0, 10)
for i in range(100000):
    h1.Fill(gRandom.Uniform(0, 10))
    h2.Fill(gRandom.Uniform(0, 10))

    if i%20 == 0:
        h1.Fill(gRandom.Gaus(3, 1))

h1.Write()
h2.Write()

oFile.Close()
