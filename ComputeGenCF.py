import yaml
import sys

from fempy import CorrelationFunction

from ROOT import TFile, gROOT


def IsLamParMatValid(lam_par):
    total = 0
    for _, h_lam in lam_par.items():
        for _, l_lam in h_lam.items():
            total += l_lam
    return abs(total - 1) < 1e-6


def SumLamPar(lam_par, treamtments):
    lam_par_summed = {
        'gen': 0,
        'flat': 0,
    }
    for hk, h_lam in lam_par.items():
        for lk, l_lam in h_lam.items():
            lam_par_summed[treamtments[hk][lk]] += l_lam
    return lam_par_summed


if __name__ == '__main__':
    inFile = TFile("~/an/DstarPi/18_fixmix/cf/RawCF_data_nopc_kStarBW15MeV_gaus_charmMassBW0.2MeV_lowKStarCharmMassBW0.4MeV_until50MeV.root")

    hCFRaw = inFile.Get("sc/hCFPurityRewFromDataMinusBkg")

    with open("/home/daniel/phsw/fempy/cfg_gencf_DstarPi.yml", "r") as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    lam_par_matr = {}
    for h_contrib in cfg['heavy']:
        h_key = list(h_contrib.keys())[0]
        h_purity = h_contrib[h_key]['purity']
        h_frac = h_contrib[h_key]['frac']
        lam_par_matr[h_key] = {}
        for l_contrib in cfg['light']:
            l_key = list(l_contrib.keys())[0]
            l_purity = l_contrib[l_key]['purity']
            l_frac = l_contrib[l_key]['frac']

            lam_par_matr[h_key][l_key] = l_frac * l_purity * h_frac * h_purity

    print('\n\n\n---> ', lam_par_matr, 'sum=', IsLamParMatValid(lam_par_matr))

    lam_par_sum = SumLamPar(lam_par_matr, cfg['treatment'])

    CFRaw = CorrelationFunction(cf=hCFRaw)
    CFGen = (CFRaw - lam_par_sum['flat'])/lam_par_sum['gen']

    gROOT.SetBatch(False)
    CFGen.get_cf().Draw()
    hCFRaw.SetLineColor(3)
    hCFRaw.Draw('same')
