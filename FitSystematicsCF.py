'''
Script to launch the fit on the correlation functions
obtained with systematic variations.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 FitSystematics.py cfg.yml

'''

import os
import argparse
import yaml
import ctypes

# from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D, TSpline3
from ROOT import TFile, gInterpreter, TH1, TH1D

from fempy import logger as log
from fempy.utils.io import Load

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--systvar', default='0')
parser.add_argument('--debug', default=False, action='store_true')
parser.add_argument('--debugfit', default=False, action='store_true')
parser.add_argument('--debugdraw', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)
if args.debugfit:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_FIT 1')
if args.debugdraw:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL_DRAW 1')

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Define the output file
oFileName = cfg['ofilename']
if cfg['suffix']:
    oFileName += f'_{cfg["suffix"]}'
oFileName += '.root'

excludeSystVars = cfg.get('exclsystvars', [])

for iSystVar in range(1, cfg['maxsystvar']+1):
    if iSystVar not in excludeSystVars:
        os.system(f"python3 fempy/FitCF.py {cfg['fitcfg']} --systvar {iSystVar}")
