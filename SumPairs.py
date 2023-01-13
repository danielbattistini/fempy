import argparse

from ROOT import TFile

parser = argparse.ArgumentParser()
parser.add_argument('inFileName')
args = parser.parse_args()

inFileName = args.inFileName
inFile = TFile(inFileName, 'update')
inFile.ls()

summedCombs = {
    'sc': ['pp', 'mm'],
    'oc': ['pm', 'mp'],
}

histClasses = [f'TH{n}{t}' for n in [1, 2, 3] for t in ['I', 'F', 'D']]


def GetHistNamesInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetObjsInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']


def GetSubdirsInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() == 'TDirectoryFile']


def GetKeyNamesInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys())]


regions = GetSubdirsInDir(inFile.Get('pp/SE'))


def SumHists(hists):
    hSum = hists[0].Clone()
    for hist in hists[1:]:
        hSum.Add(hist)
    return hSum


for combName, combsToSum in summedCombs.items():
    for event in ['SE', 'ME']:
        inFile.mkdir(f'{combName}/{event}')
        inFile.cd(f'{combName}/{event}')

        for histoName in GetHistNamesInDir(inFile.Get('pp/SE')):
            SumHists([inFile.Get(f'{comb}/{event}/{histoName}') for comb in combsToSum]).Write()

        for region in regions:
            inFile.mkdir(f'{combName}/{event}/{region}')
            inFile.cd(f'{combName}/{event}/{region}')

            for histoName in GetHistNamesInDir(inFile.Get('pp/SE/sgn')):
                SumHists([inFile.Get(f'{comb}/{event}/{region}/{histoName}') for comb in combsToSum]).Write()
inFile.Close()
