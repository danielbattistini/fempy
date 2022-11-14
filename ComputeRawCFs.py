import uproot
import numpy as np
import sys
import argparse
import yaml
from rich import print

from ROOT import TH1D, TH2F

from fempy.utils.format import TranslateToLatex
from fempy.utils.analysis import is_mass_selected, Pair

from fempy import CorrelationFunction

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--pair')
parser.add_argument('in_file')
parser.add_argument('o_file')
args = parser.parse_args()

in_file_name = args.in_file
o_file_name = args.o_file

pair = Pair(args.pair)

bin_widths = [1, 2, 4, 5, 10, 20, 40, 50]  # MeV/c


def get_list_of_dirs(path):
    '''Get the list of the syst variations in the input file'''
    return [key[:-2] for key in path.keys() if '/' not in key]


with uproot.recreate(o_file_name) as o_file:
    hist = np.histogram(np.random.normal(0, 1, 1))
    hist2 = np.histogram(np.random.normal(0, 1, 1), np.random.normal(0, 1, 1))
    o_file['hDummy'] = hist
    o_file['hDummy2'] = hist2
    del o_file['hDummy']
    del o_file['hDummy2']

    with uproot.open(in_file_name) as in_file:
        d_se, d_me, d_cf = ({} for _ in range(3))

        variations = get_list_of_dirs(in_file['distr'])
        print(variations)
        variations = [
            '0',
            # 'oldpckept',
            # 'oldpckept_motherpi_eq413',
            # 'newpckept_motherpi_eq413',
            # 'oldpckept_motherpi_neq413',
            # 'newpckept_motherpi_neq413',
            # 'oldpcrm',
            # 'oldpcrm_motherpi_eq413',
            # 'newpcrm_motherpi_eq413',
            # 'oldpcrm_motherpi_neq413',
            # 'newpcrm_motherpi_neq413',
            # 'oldpcrm_multDeq1',
            # 'oldpcrm_multDgt1',
            # 'newpcrm_multDeq1',
            # 'newpcrm_multDgt1'
        ]
        # load SE and ME
        for mass_region in pair.mass_regions:
            for var in variations:
                for comb in ['pp', 'mm', 'pm', 'mp']:
                    d_se[f'distr/{var}/hSE_{comb}_{mass_region}'] = in_file[f'distr/{var}/hSE_{comb}_{mass_region}'].to_pyroot()
                    d_me[f'distr/{var}/hME_{comb}_{mass_region}'] = in_file[f'distr/{var}/hME_{comb}_{mass_region}'].to_pyroot()

                # sum SE same charge
                distr_pp = d_se[f'distr/{var}/hSE_pp_{mass_region}'].Clone()
                distr_mm = d_se[f'distr/{var}/hSE_mm_{mass_region}']
                distr_sc = distr_pp.Clone()
                distr_sc.Add(distr_mm)
                distr_sc.SetDirectory(0)
                d_se[f'distr/{var}/hSE_sc_{mass_region}'] = distr_sc

                # sum SE opposite charge
                distr_pm = d_se[f'distr/{var}/hSE_pm_{mass_region}'].Clone()
                distr_mp = d_se[f'distr/{var}/hSE_mp_{mass_region}']
                distr_oc = distr_pm.Clone()
                distr_oc.Add(distr_mp)
                distr_oc.SetDirectory(0)
                d_se[f'distr/{var}/hSE_oc_{mass_region}'] = distr_oc

                # sum ME same charge
                distr_pp = d_me[f'distr/{var}/hME_{comb}_{mass_region}'].Clone()
                distr_mm = d_me[f'distr/{var}/hME_{comb}_{mass_region}']
                distr_sc = distr_pp.Clone()
                distr_sc.Add(distr_mm)
                distr_sc.SetDirectory(0)
                d_me[f'distr/{var}/hME_sc_{mass_region}'] = distr_sc

                # sum ME opposite charge
                distr_pm = d_me[f'distr/{var}/hME_{comb}_{mass_region}'].Clone()
                distr_mp = d_me[f'distr/{var}/hME_{comb}_{mass_region}']
                distr_oc = distr_pm.Clone()
                distr_oc.Add(distr_mp)
                distr_oc.SetDirectory(0)
                d_me[f'distr/{var}/hME_oc_{mass_region}'] = distr_oc

        # project SE
        d_se_proj = {}
        for distr in d_se:
            proj = d_se[distr].ProjectionX(f'{distr}', 0, d_se[distr].GetNbinsX())
            d_se_proj[f'{distr}'] = proj

        # project ME
        d_me_proj = {}
        for distr in d_me:
            proj = d_me[distr].ProjectionX(f'{distr}_unrew', 0, d_me[distr].GetNbinsX())
            d_me_proj[f'{distr}_unrew'] = proj

        for mass_region in pair.mass_regions:
            for comb in ['pp', 'mm', 'sc', 'pm', 'mp', 'oc']:
                for bw in bin_widths:
                    for var in variations:
                        se = d_se_proj[f'distr/{var}/hSE_{comb}_{mass_region}'].Clone()
                        me = d_me_proj[f'distr/{var}/hME_{comb}_{mass_region}_unrew'].Clone()

                        se.Rebin(bw)
                        me.Rebin(bw)

                        d_cf[f'distr/{var}/hSE_{comb}_{mass_region}_bw{bw}MeV'] = se
                        d_cf[f'distr/{var}/hME_{comb}_{mass_region}_unrew_bw{bw}MeV'] = me

                        cf = CorrelationFunction(se, me, pair.norm_range, um='MeV').get_cf()
                        d_cf[f'cf/{var}/hCF_{comb}_{mass_region}_unrew_bw{bw}MeV'] = cf

        for key, distr in {**d_se, **d_me}.items():
            o_file[key] = distr

        for cf_key in d_cf:
            o_file[cf_key] = d_cf[cf_key]

        # compute multiplicity reweighted CFs
        mult_bins_mins = list(range(0, 180, 5))
        mult_bins_maxs = list(range(5, 185, 5))

        for mass_region in pair.mass_regions:
            for comb in ['pp', 'mm', 'sc', 'pm', 'mp', 'oc']:
                for var in variations:
                    print(f'distr/{var}/weights/hSE_{comb}_{mass_region}')
                    se = d_se[f'distr/{var}/hSE_{comb}_{mass_region}']
                    print(se, "\n\n")
                    mult_se = se.ProjectionY(f'distr/{var}/hSE_{comb}_{mass_region}_mult', 0, se.GetNbinsY()+1)
                    mult_se.Rebin(5)
                    o_file[f'distr/{var}/weights/hSE_{comb}_{mass_region}_mult'] = mult_se

                    me = d_me[f'distr/{var}/hME_{comb}_{mass_region}']
                    mult_me = me.ProjectionY(f'distr/{var}/weights/hME_{comb}_{mass_region}_mult', 0, me.GetNbinsY()+1)
                    mult_me.Rebin(5)
                    o_file[f'distr/{var}/weights/hME_{comb}_{mass_region}_mult'] = mult_me

                    # compute weights
                    h_weight = mult_se.Clone('hMultWeights')
                    weights = []
                    me_mult = []
                    for i_bin in range(0, h_weight.GetNbinsX()+1):
                        me_counts = mult_me.GetBinContent(i_bin)
                        if me_counts == 0:
                            weight = 0
                        else:
                            weight = mult_se.GetBinContent(i_bin+1)/me_counts

                        weights.append(weight)
                        h_weight.SetBinContent(i_bin, weight)

                    o_file[f'cf/{var}/weights/hME_{comb}_{mass_region}_multweights'] = h_weight

                    # compute mult slices
                    me_proj = []
                    for (mult_min, mult_max) in zip(mult_bins_mins, mult_bins_maxs):
                        first_mult_bin = me.GetXaxis().FindBin(mult_min + 1.e-6)
                        last_mult_bin = me.GetXaxis().FindBin(mult_max - 1.e-6)

                        me_proj.append(me.ProjectionX(f'{distr}_mult{mult_min}_{mult_max}', first_mult_bin, last_mult_bin))

                        # o_file[f'cf/{var}/hME_{comb}_{mass_region}_mult{mult_min}_{mult_max}'] = proj

                    # compute reweighted ME
                    h_me_rew = TH1D("", "", me_proj[0].GetNbinsX(), 0, me_proj[0].GetXaxis().GetXmax())
                    for (me, w) in zip(me_proj, weights):
                        h_me_rew.Add(me, w)
                    for bw in bin_widths:
                        h_me_rew_rebin = h_me_rew.Clone()
                        h_me_rew_rebin.Rebin(bw)
                        o_file[f'distr/{var}/hME_{comb}_{mass_region}_rew_bw{bw}MeV'] = h_me_rew_rebin

                        # compute the reweighted CF
                        se = d_cf[f'distr/{var}/hSE_{comb}_{mass_region}_bw{bw}MeV']

                        cf = CorrelationFunction(se, h_me_rew_rebin, pair.norm_range, um='MeV').get_cf()
                        o_file[f'cf/{var}/hCF_{comb}_{mass_region}_rew_bw{bw}MeV'] = cf
print(f"Correlation functions saved in {o_file_name}")
