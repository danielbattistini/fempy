from ROOT import TFile, TCanvas, gInterpreter, TF1, TDatabasePDG
gInterpreter.ProcessLine('#include "CorrelationFitter.hxx"')
from ROOT import CorrelationFitter, BreitWigner

from fempy.utils.io import Load

for eta in ['mid', 'fwd']:
    for level in ['gen', 'reco']:
        inFile = TFile(f'/home/ktas/ge86rim/an/alice3/femto/cf/RawCF_KK_{level}_B-1.0T_geom-20230731_pileup-1_hf_forced-decays_{eta}_cpr.root')
        oFile = TFile(f'/home/ktas/ge86rim/an/alice3/femto/cf/FitRawCF_KK_{level}_B-1.0T_geom-20230731_pileup-1_hf_forced-decays_{eta}_cpr.root', 'recreate')

        hCF = Load(inFile, 'p03_12/sgn/hCF')
        fitter = CorrelationFitter(hCF, 0, 0.3)
        fitter.Add('pol2', [('p0', 1, 0.5, 1.5), ('p1', 0, -0.5, 0.5), ('p2', 0, -1, 5)])
        fitter.Add('bw', [
            ('yield', 0.03, 0.001, 0.5),
            ('mean', 0.126, 0.125, 0.127),
            ('gamma', 0.004, 0.003, 0.05),
            ])
        fitter.Fit()
        cFit = TCanvas('cFit', '', 600, 600)
        fitter.Draw(cFit)
        cFit.SaveAs(f'/home/ktas/ge86rim/an/alice3/femto/cf/FitRawCF_KK_{level}_B-1.0T_geom-20230731_pileup-1_hf_forced-decays_{eta}_cpr.pdf')

        hCF.Write()
        fitter.GetFunction().Write()
