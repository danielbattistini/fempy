import argparse

from ROOT import TFile

from fempy.utils.io import GetPairName, GetObjectFromFile

def CountPairs(inFileName, fd=False):
    inFile = TFile.Open(inFileName)

    if fd:
        inFile.ls()
        names = GetPairName(inFile)
        for name in names:
            print(name)
            FDhistBasePath = f"HM_CharmFemto_{name}_Results0/HM_CharmFemto_{name}_Results0"
            
            # same charge
            nEntries = 0
            nEntriesLowKStar = 0
            for p1, p2 in [(0, 2), (1, 3)]:

                hist = GetObjectFromFile(inFile, f"{FDhistBasePath}/Particle{p1}_Particle{p2}/SEDist_Particle{p1}_Particle{p2}")
                nEntries += hist.GetEntries()
                nEntriesLowKStar += hist.Integral(1, hist.GetXaxis().FindBin(0.1999999))
            print('nEntries in sc: Full range:', nEntries, 'low k*: ', nEntriesLowKStar)
            
            # opposite charge
            nEntries = 0
            nEntriesLowKStar = 0
            for p1, p2 in [(0, 3), (1, 2)]:
                hist = GetObjectFromFile(inFile, f"{FDhistBasePath}/Particle{p1}_Particle{p2}/SEDist_Particle{p1}_Particle{p2}")
                nEntries += hist.GetEntries()
                nEntriesLowKStar += hist.Integral(1, hist.GetXaxis().FindBin(0.1999999))
            print('nEntries in oc: Full range:', nEntries, 'low k*: ', nEntriesLowKStar)

    else:
        for comb in ['pp', 'mm', 'pm', 'mp', 'sc', 'oc']:
            hKStarVsMult = inFile.Get(f'{comb}/SE/sgn/hMultVsKStar0')
            hKStar = hKStarVsMult.ProjectionX()
            hKStar.Draw()
            nEntries = hKStar.GetEntries()
            nEntriesLowKStar = hKStar.Integral(1, hKStar.GetXaxis().FindBin(199.9999))
            print(f'nEntries in {comb}: Full range:', nEntries, 'low k*: ', nEntriesLowKStar)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument('--fd', action='store_true', help='FemtoDream', default=False)
    args = parser.parse_args()

    CountPairs(args.inFile, args.fd)
