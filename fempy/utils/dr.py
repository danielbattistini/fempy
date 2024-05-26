#!/usr/bin/env python
'''
Download AnalysisResults.root from the grid

Install running install.sh
'''

import os
import subprocess
import argparse
import sys

from fempy import logger as log

parser = argparse.ArgumentParser()
parser.add_argument('train', help='Name of the train from which file should be downloaded')
parser.add_argument('run', help='Number of the run')
args = parser.parse_args()

train = args.train
run = args.run

# Check that the train is valid
if 'CF' in train:
    pwg = 'CF'
elif 'D2H' in train:
    pwg = 'HF'
elif 'NanoAOD' in train:
    pwg = 'ZZ'
else:
    log.critical('Unknown train')

# Search the files
path = f'/alice/cern.ch/user/a/alitrain/PWG{pwg}/{train}/{run}_.*'
try:
    files = subprocess.check_output(f'alien_find -r "{path}/AnalysisResults.root"', shell=True, universal_newlines='\n')
except subprocess.CalledProcessError:
    log.critical('No files found in %s', path)
files = str(files).split('\n')
files.remove('')

children = []
merged = []
for name in files:
    if 'child' in name:
        children.append(name)
    else:
        merged.append(name)

print('Found the following files:\nmerged:')
for name in merged:
    print('   ', name)
print('\nChildren:')
for name in children:
    print('   ', name)

choice = None
while choice not in ['m', 'c']:
    choice = input('Download type: merged (m) or children (c): ')

if choice == 'm':
    if len(merged) > 1:
        sys.exit()
    os.system(f'alien_cp {merged[0]} file:AnalysisResults_{run}.root')

if choice == 'c':
    for name in children:
        start = name.index('child_')+6
        child = name[start:start+name[start:].index('/')]

        os.system(f'alien_cp {name} file:merge_{run}/child_{child}/AnalysisResults.root')
