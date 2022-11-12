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
        print(variations)
        variations = [
                    '0',
            'oldpckept',
            'oldpckept_motherpi_eq413',
            'newpckept_motherpi_eq413',
            'oldpckept_motherpi_neq413',
            'newpckept_motherpi_neq413',
            'oldpcrm',
            'oldpcrm_motherpi_eq413',
            'newpcrm_motherpi_eq413',
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
            proj = d_se[distr].ProjectionX(f'{distr}_unrew', 0, d_se[distr].GetNbinsX())
            d_se_proj[f'{distr}_unrew'] = proj

        # project ME
        d_me_proj = {}
        for distr in d_me:
            proj = d_me[distr].ProjectionX(f'{distr}_unrew', 0, d_me[distr].GetNbinsX())
            d_me_proj[f'{distr}_unrew'] = proj

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
                    for i_bin in range(0, h_weight.GetNbinsX()+1):
                        me_counts = mult_me.GetBinContent(i_bin)
                        if me_counts == 0:
                            h_weight.SetBinContent(i_bin, 0)
                        else:
                            h_weight.SetBinContent(i_bin, mult_se.GetBinContent(i_bin+1)/me_counts)
                    o_file[f'distr/{var}/weights/hME_{comb}_{mass_region}_multweights'] = h_weight


        # hWeight.Divide(mult_me)
        # make mult slices

        # d_me_proj = {}
        # for distr in d_me:
        #     for (mult_min, mult_max) in zip(mult_bins_mins, mult_bins_maxs):
        #         first_mult_bin = d_me[distr].GetXaxis().FindBin(mult_min)
        #         last_mult_bin = d_me[distr].GetXaxis().FindBin(mult_max)

        #         proj = d_me[distr].ProjectionX(f'{distr}_mult{mult_min}_{mult_max}', first_mult_bin, last_mult_bin)
        #         d_me_proj[f'{distr}_mult{mult_min}_{mult_max}'] = proj

        # print(me)
        # # make mult projection
        # mult_bins_mins = list(range(0, 180, 5))
        # mult_bins_maxs = list(range(5, 185, 5))

        # # project SE
        # d_se_proj = {}
        # for distr in d_se:
        #     print(d_se[distr])
        #     for (mult_min, mult_max) in zip(mult_bins_mins, mult_bins_maxs):
        #         first_mult_bin = d_se[distr].GetXaxis().FindBin(mult_min)
        #         last_mult_bin = d_se[distr].GetXaxis().FindBin(mult_max)

        #         proj = d_se[distr].ProjectionX(f'{distr}_mult{mult_min}_{mult_max}', first_mult_bin, last_mult_bin)
        #         d_se_proj[f'{distr}_mult{mult_min}_{mult_max}'] = proj

        # # project ME
        # d_me_proj = {}
        # for distr in d_me:
        #     print(d_me[distr])
        #     for (mult_min, mult_max) in zip(mult_bins_mins, mult_bins_maxs):
        #         first_mult_bin = d_me[distr].GetXaxis().FindBin(mult_min)
        #         last_mult_bin = d_me[distr].GetXaxis().FindBin(mult_max)

        #         proj = d_me[distr].ProjectionX(f'{distr}_mult{mult_min}_{mult_max}', first_mult_bin, last_mult_bin)
        #         d_me_proj[f'{distr}_mult{mult_min}_{mult_max}'] = proj

        # def is_check(name):
        #     for check in checks:
        #         if check in name:
        #             return True, check
        #     return False, None

        # for d_distr in[d_se, d_me, d_se_proj, d_me_proj, d_cf]:

        # for d_distr in [d_se, d_me, d_se_proj, d_me_proj]:
        #     for distr_key, distr in d_distr.items():
        #         (is_a_check, check_name) = is_check(distr_key)
        #         if is_a_check:
        #             print(distr_key)
        #             o_file[f'checks/distr/{check_name}/{distr_key[:2]}'] = distr
        #         else:
        #             o_file[f'distr/{distr_key}'] = distr

        # for distr_key, distr in d_cf.items():
        #     (is_a_check, check_name) = is_check(distr_key)
        #     if is_a_check:
        #         o_file[f'checks/cf/{check_name}/{distr_key}'] = distr
        #     else:
        #         o_file[f'cf/{distr_key}'] = distr

            # mult compute projections
            # distr.ProjectionX()
            # nBins = round(pair.max_kstar / bw)

            # cf = CorrelationFunction(se, me)
            #                 cf = CorrelationFunction(se, me, pair.norm_range)
            #                 in_file[f'bw{bw}MeV/hCF_{comb}_{mass_region}'] = cf.get_cf()

            #         # perform checks
            #         for (check_name, check_sel) in pair.checks:

            #             for bw in bin_widths:
            #                 nBins = round(pair.max_kstar / bw)
            #                 h_kstar = TH1D(f'h{event}_{comb}_{mass_region}_{bw}_{check_name}', "", nBins, 0, pair.max_kstar)
            #                 print(f'h{event}_{comb}_{mass_region}_{bw}_{check_name}', h_kstar)

            #                 in_file[f'checks/{check_name}/bw{bw}MeV/h{event}_{comb}_{mass_region}'] = h_kstar

            #         # quality assurance:: mass regions
            #         if 'Dstar' in pair.name:
            #             mass_lower = 0.12
            #             mass_upper = 0.25

            #         h_mass_pthf = TH2F(f"hMass_PtHF_{comb}_{mass_region}", "", 200, 0, 10, 200, mass_lower, mass_upper)
            #         h_mass_pthf.FillN(len(pt), pt, mass, np.ones_like(pt, 'd'))
            #         in_file[f'qa/{comb}_{mass_region}/h{event}_Mass_PtHF'] = h_mass_pthf

            #         # multiplicity
            #         mult_hf = np.array(df_proj['heavy_mult'], 'd')
            #         mult_lf = np.array(df_proj['light_mult'], 'd')
            #         h_multhf_multlf = TH2F(f"hMultHF_MultLF_{comb}_{mass_region}", "", 16, -0.5, 15.5, 31, -0.5, 30.5)
            #         h_multhf_multlf.FillN(len(mult_hf), mult_hf, mult_lf, np.ones_like(mult_hf, 'd'))
            #         in_file[f'qa/{comb}_{mass_region}/h{event}_hMultHF_MultLF'] = h_multhf_multlf

            #     # sum charge combinaitons
            #     combs_to_sum = [
            #         ('sc', ['pp', 'mm']),
            #         ('oc', ['pm', 'mp'])
            #     ]
            #     for mass_region in pair.mass_regions:
            #         for bw in bin_widths:
            #             for comb_to_sum in combs_to_sum:
            #                 (derived_comb, primitive_combs) = comb_to_sum

            #                 # central selection (0) and syst variation
            #                 for i_syst, syst_sel, in enumerate(pair.syst_sel):
            #                     sum_que = []
            #                     for comb in primitive_combs:
            #                         sum_que.append(in_file[f'bw{bw}MeV/h{event}_{comb}_{mass_region}'].to_hist())
            #                     h_sum = sum_que[0]
            #                     for hist in sum_que[1:]:
            #                         h_sum += hist
            #                     in_file[f'bw{bw}MeV/h{event}_{derived_comb}_{mass_region}'] = h_sum

            #                 # checks
            #                 for (check_name, check_sel) in pair.checks:
            #                     sum_que = []
            #                     for comb in primitive_combs:
            #                         sum_que.append(in_file[f'checks/{check_name}/bw{bw}MeV/h{event}_{comb}_{mass_region}'].to_hist())
            #                     h_sum = sum_que[0]
            #                     for hist in sum_que[1:]:
            #                         h_sum += hist
            #                     in_file[f'checks/{check_name}/bw{bw}MeV/h{event}_{derived_comb}_{mass_region}'] = h_sum

            # # compute the correlation functions
            # for mass_region in pair.mass_regions:
            #     for bw in bin_widths:
            #         for comb in ['pp', 'mm', 'sc', 'pm', 'mp', 'oc']:
            #             for i_syst, _, in enumerate(pair.syst_sel):
            #                 se = in_file[f'bw{bw}MeV/hSE_{comb}_{mass_region}']
            #                 me = in_file[f'bw{bw}MeV/hME_{comb}_{mass_region}']
            #                 cf = CorrelationFunction(se, me, pair.norm_range, um='MeV').get_cf()
            #                 in_file[f'bw{bw}MeV/hCF_{comb}_{mass_region}'] = cf

            #             for (check_name, _) in pair.checks:
            #                 se = in_file[f'checks/{check_name}/bw{bw}MeV/hSE_{comb}_{mass_region}']
            #                 me = in_file[f'checks/{check_name}/bw{bw}MeV/hME_{comb}_{mass_region}']
            #                 cf = CorrelationFunction(se, me, pair.norm_range, um='MeV').get_cf()
            #                 in_file[f'checks/{check_name}/bw{bw}MeV/hCF_{comb}_{mass_region}'] = cf
print(f"Correlation functions saved in {o_file_name}")
