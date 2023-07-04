#!/home/daniel/alice/sw/ubuntu2004_x86-64/Python/v3.9.12-4/bin/python3
'''
Computes the mean-proper decay lenght (ctau) in fm of a particle given the width in MeV
'''

import argparse

hbarc = 197.3269804 # MeV * fm

parser = argparse.ArgumentParser()
parser.add_argument('width', type=float, help='(MeV)')
args = parser.parse_args()

width = args.width
print(f'ctau = {hbarc / width:.2f} fm')