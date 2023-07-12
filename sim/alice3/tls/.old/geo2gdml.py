from math import *
import ROOT
ROOT.TGeoManager.Import("o2sim_geometry.root")
ROOT.gGeoManager.Export("o2sim_geometry.gdml")
