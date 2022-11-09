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