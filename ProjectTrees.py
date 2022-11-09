import uproot
import numpy as np
import argparse

from ROOT import TH1D, TH2F

from fempy.utils.analysis import is_mass_selected, Pair

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('in_file')
parser.add_argument('o_file')
parser.add_argument('--pair')
parser.add_argument('--redstat', action='store_true', default=False)
args = parser.parse_args()

in_file_name = args.in_file
o_file_name = args.o_file

pair = Pair(args.pair)

bin_widths = [1, 2, 4, 5, 10, 20, 40, 50]  # MeV/c
cols_to_keep = [
    'kStar',
    'mult',
    'heavy_mult',
    'heavy_invmass',
    'heavy_pt',
    'is_newpcrm',
    'is_oldpcrm',
    'light_eta',
    'light_px',
    'light_py',
    'light_mult',
    'light_nsigtpc',
    'light_nsigtof',
    'light_dcaxy',
    'light_dcaz',
    'light_ncrossed',
    'light_ncls',
    'light_motherpdg',
]

with uproot.recreate(o_file_name) as o_file:
    # hist = np.histogram(np.random.normal(0, 1, 1))
    # o_file['hDummy'] = hist
    # del o_file['hDummy']
    # o_file['qa/hist'] = hist
    # del o_file['qa/hist']

    print(f'Opening {in_file_name} ... ', end='')
    with uproot.open(f'{in_file_name}:HM_CharmFemto_DstarPion_Trees0') as dir:
        print('done!')

        for event in ['SE', 'ME']:
            for comb in ['pp', 'mm', 'pm', 'mp']:
                tree_name = f't{event}_{comb}'
                print(f'\nloading {tree_name} ... ', end='')
                df = dir[tree_name].arrays(library='pd', filter_name=cols_to_keep)
                print('done!')

                if args.redstat:
                    df = df.head(10000)

                df.loc[:, 'kStar'] *= 1000  # convert GeV/c to MeV/c

                for mass_region in pair.mass_regions:
                    print(f'    Mass region: {mass_region}')
                    df_proj = df[df.apply(lambda row: is_mass_selected(row['heavy_invmass'], row['heavy_pt'],
                                                                       pair.hpdg, mass_region), axis=1)]

                    hist_basename = f'h{event}_{pair.name}_{comb}_{mass_region}'
                    for i_syst, syst_sel, in enumerate(pair.syst_sel):
                        print(f'        syst var: {i_syst}')

                        df_syst = df_proj.query(syst_sel)
                        values = np.array(df_syst['kStar'], 'd')

                        if len(values) == 0:
                            print(f"\033[33m *** {hist_basename}_{i_syst} has 0 entries!\033[0m")
                        for bw in bin_widths:
                            nBins = round(pair.max_kstar / bw)
                            h_kstar = TH1D(f'{hist_basename}_{bw}_{i_syst}', "", nBins, 0, pair.max_kstar)

                            if len(values) > 0:
                                h_kstar.FillN(len(values), values, np.ones_like(values, 'd'))
                            o_file[f'bw{bw}MeV/{hist_basename}_{i_syst}'] = h_kstar

                    # perform checks
                    for (check_name, check_sel) in pair.checks:
                        print(f'        check: {check_name}')
                        df_check_sel = df_proj.query(check_sel)
                        values = np.array(df_check_sel['kStar'], 'd')
                        if len(values) == 0:
                            print(f"\033[33m *** {check_name} {hist_basename} has 0 entries!\033[0m")
                        for bw in bin_widths:
                            nBins = round(pair.max_kstar / bw)
                            h_kstar = TH1D(f'{hist_basename}_{bw}_{check_name}', "", nBins, 0, pair.max_kstar)
                            # print(f'{hist_basename}_{bw}_{check_name}', h_kstar)

                            if len(values) > 0:
                                h_kstar.FillN(len(values), values, np.ones_like(values, 'd'))
                            o_file[f'checks/{check_name}/bw{bw}MeV/{hist_basename}'] = h_kstar

                    # quality assurance: mass regions
                    if 'Dstar' in pair.name:
                        mass_lower = 0.12
                        mass_upper = 0.25
                    pt = np.array(df_proj['heavy_pt'], 'd')
                    mass = np.array(df_proj['heavy_invmass'], 'd')

                    h_mass_pthf = TH2F(f"hMass_PtHF_{comb}_{mass_region}", "", 200, 0, 10, 200, mass_lower, mass_upper)
                    h_mass_pthf.FillN(len(pt), pt, mass, np.ones_like(pt, 'd'))
                    o_file[f'qa/{comb}_{mass_region}/h{event}_Mass_PtHF'] = h_mass_pthf

                    # multiplicity
                    mult_hf = np.array(df_proj['heavy_mult'], 'd')
                    mult_lf = np.array(df_proj['light_mult'], 'd')
                    h_multhf_multlf = TH2F(f"hMultHF_MultLF_{comb}_{mass_region}", "", 16, -0.5, 15.5, 31, -0.5, 30.5)
                    h_multhf_multlf.FillN(len(mult_hf), mult_hf, mult_lf, np.ones_like(mult_hf, 'd'))
                    o_file[f'qa/{comb}_{mass_region}/h{event}_hMultHF_MultLF'] = h_multhf_multlf

                    # mult vs kstar
                    mult = np.array(df_proj['mult'], 'd')
                    k_star = np.array(df_proj['kStar'], 'd')
                    h_kstar_mult = TH2F(f"hKStar_Mult_{comb}_{mass_region}", "", pair.max_kstar, 0, pair.max_kstar, 180, 0.5, 180.5)
                    h_kstar_mult.FillN(len(mult), k_star, mult, np.ones_like(k_star, 'd'))
                    o_file[f'qa/{comb}_{mass_region}/h{event}_hKStar_Mult'] = h_kstar_mult
