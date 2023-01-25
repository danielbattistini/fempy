import os
import sys

from ROOT import TFile, TDirectoryFile, TList

def GetObjectFromFile(inFile, pathToObj):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TFile
     - TDirecotryFile
     - TList

    Parameters
    -----------
    inFile: TFile of the input file
    pathToObj: path of the object inside the root file

    Returns:
    -----------
    outObj: target root object
    '''

    pathToObj = os.path.normpath(pathToObj)
    pathElements = pathToObj.split(os.sep)
    outObj = inFile.Get(pathElements.pop(0))

    for _, containerName in enumerate(pathElements):
        if isinstance(outObj, (TFile, TDirectoryFile)):
            outObj = outObj.Get(containerName)
        elif isinstance(outObj, TList):
            outObj = outObj.FindObject(containerName)
        else:
            print(f'\033[31mError\033[0m: instance of {type(outObj)} not implemented. Exit!')
            sys.exit()

    return outObj


histClasses = [f'TH{n}{t}' for n in [1, 2, 3] for t in ['I', 'F', 'D']]
def GetHistNamesInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetObjsInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']


def GetSubdirsInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys()) if key.GetClassName() == 'TDirectoryFile']


def GetKeyNamesInDir(directory):
    return [key.GetName() for key in list(directory.GetListOfKeys())]
