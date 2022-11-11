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
    o_file['hDummy'] = hist
    del o_file['hDummy']
    with uproot.open(in_file_name) as in_file:
        d_se, d_me, d_cf = ({} for _ in range(3))

        variations = get_list_of_dirs(in_file['distr'])

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
            proj = d_se[distr].ProjectionX(f'{distr}_unrew', 0, d_se[distr].GetNbinsX())
            d_se_proj[f'{distr}_unrew'] = proj

        # project ME
        d_me_proj = {}
        for distr in d_me:
            proj = d_me[distr].ProjectionX(f'{distr}_unrew', 0, d_me[distr].GetNbinsX())
            d_me_proj[f'{distr}_unrew'] = proj

        print(d_se)
        print(d_me)
        for mass_region in pair.mass_regions:
            for comb in ['pp', 'mm', 'sc', 'pm', 'mp', 'oc']:
                for bw in bin_widths:
                    for var in variations:
                        se = d_se_proj[f'distr/{var}/hSE_{comb}_{mass_region}_unrew'].Clone()
                        me = d_me_proj[f'distr/{var}/hME_{comb}_{mass_region}_unrew'].Clone()

                        se.Rebin(bw)
                        me.Rebin(bw)

                        d_cf[f'distr/{var}/hSE_{comb}_{mass_region}_unrew_bw{bw}MeV'] = se
                        d_cf[f'distr/{var}/hME_{comb}_{mass_region}_unrew_bw{bw}MeV'] = me

                        cf = CorrelationFunction(se, me, pair.norm_range, um='MeV').get_cf()
                        d_cf[f'cf/{var}/hCF_{comb}_{mass_region}_unrew_bw{bw}MeV'] = cf

        for key, distr in {**d_se, **d_me}.items():
            print(key)
            o_file[key] = distr

        for cf_key in d_cf:
            print(cf_key)
            o_file[cf_key] = d_cf[cf_key]

print(f"Correlation functions saved in {o_file_name}")
