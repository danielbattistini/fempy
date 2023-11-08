#!/bin/bash
# Script to install fempy

mkdir -p ~/.local/bin/

cp fempy/utils/width2ctau.py ~/.local/bin/width2ctau
chmod +x ~/.local/bin/width2ctau

cp fempy/utils/dr.py ~/.local/bin/dr
chmod +x ~/.local/bin/dr

pip3 install .
